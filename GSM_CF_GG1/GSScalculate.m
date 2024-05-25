% -------------------------------------------------------------------------
% Функция расчета
% -------------------------------------------------------------------------

function [lambda, j, t, to, V, f, No, ca, ro, Bl] = GSScalculate(M, N, Q, MU, CS, W, errorRate)

    W = getW(M, Q); % матрица вероятностей нахождения заявки в текущем узле

    mostLoadedNode = getMostLoadedNode(M, MU, W); % наиболее загруженный узел

    l = 0;
    x = [];
    y = [];
    
    Bl = (MU(mostLoadedNode) / W(mostLoadedNode)) * (1 - 1 / N);
    
    x = [x Bl];
    y = [y Bl];
    
    L = 0;
    
    while abs(L - N) >= errorRate
        l = l + 1;
        [lambda, j, t, to] = doIteration(Bl);
        j = transpose(j);
        L = sum(j);
        Bl = Bl * (N / L);
        y = [y Bl];
        if l > 1
            Bl = getBlWegstein(l, x, y);
        end
        x = [x Bl];
    end

    V = getV(lambda);
    f = getF(V);
    No = getNo(j, to, t);

% -------------------------------------------------------------------------
% Вспомогательные функции
% -------------------------------------------------------------------------
    
    % lambda_ij
    
    function lambda_ij = getLambdaij(lambda)
        lambda_ij = zeros(M);
        
        for i = 1:M
            for j = 1:M
                lambda_ij(i, j) = lambda(i) * Q(i, j);
            end
        end
    end
    
    % ca0
    
    function [ca] = ca0(M)
        ca = zeros(M,1);
        ca(1)=1;
    end

    % ca

    function [ca, ro] = getCai(lambda, lambda_ij)
        ro = zeros(M,1);
    
        for i = 1:M
            ro(i) = lambda(i) / MU(i);
        end
        
        function [F] = caFun(ca)
            Wi = zeros(M,1);
            Vi = zeros(M,1);
            F = zeros(M,1);
            
            for i = 1:M
                tmp = 0;
            
                for j = 1:M
                    tmp = tmp + (lambda_ij(j, i) / lambda(i))^2;
                end
            
                Vi(i) = 1 / tmp;
                Wi(i) = 1 / (1 + 4 * (1 - ro(i))^2 * (Vi(i) - 1));
            end
            
            for i = 1:M
                tmp = 0;
            
                for j = 1:M
                    ca_ij = Q(j, i) * (ro(j)^2 * CS(j) + (1 - ro(j)^2) * ca(j)) + 1 - Q(j,i);
                    tmp = tmp + lambda_ij(j, i) * ca_ij;
                end
                tmp = tmp / lambda(i);
                F(i) = Wi(i) * tmp + 1 - Wi(i) - ca(i);
            end
        
        end
        
        ca = fsolve(@caFun, ca0(M));
    end
    
    % to без корректирующего фактора
    
    function toWCF = getToWCF(ro, ca)
        toWCF = zeros(M,1);
        
        for i = 1:M
            ro_i = ro(i);
            ca_i = ca(i);
            to = ro_i * (CS(i) + ca_i) / (2 * MU(i) * (1 - ro_i));
            if ca_i < 1
                to = to * exp((ro_i - 1) * ((1 - ca_i)^2) / (1.5 * ro_i * (CS(i) + ca_i)));
            else
                to = to * exp((ro_i - 1) * (ca_i - 1) / (4 * CS(i) + ca_i));
            end
            toWCF(i) = to;
        end
    end
    
    % Одна итерация алгоритма
    
    function [lambda, j, t, to] = doIteration(Bl)
        lambda = Bl * W;
        lambda_ij = getLambdaij(lambda);
        [ca, ro] = getCai(lambda, lambda_ij);
        toWCF = getToWCF(ro, ca);
        
        to = zeros(M,1);
        t = zeros(M,1);
        j = zeros(M,1);
        
        for i = 1:M
            to_i = toWCF(i) * (N - 1) / (N + toWCF(i) * MU(i));
            to(i) = to_i;
            t(i) = to_i + 1 / MU(i);
            j(i) = t(i) * lambda(i);
        end
    end
    
    % Метод Вегстейна
    
    function Bl_W = getBlWegstein(l, x, y)
        xk = x((l - 1):l);
        yk = y(l:(l + 1));
        
        Bl_W = yk(2) - ((yk(2) - yk(1)) * (yk(2) - xk(2))) / (yk(2) - yk(1) - xk(2) + xk(1));
    end
    
    % Время цикла заявок для каждого из узлов сети
    
    function V = getV(lambda)
        V = zeros(1,M);
        
        for i = 1:M
            V(i) = N / lambda(i);
        end
    end
    
    % Пропускная способность сети
    
    function f = getF(V)
        f = M / sum(V);
    end
    
    % Число заявок в очереди
    
    function No = getNo(j, to, t)
        No = zeros(1,M);
    
        for i = 1:M
            No(i) = j(i) * to(i) / t(i);
        end
    end
end

