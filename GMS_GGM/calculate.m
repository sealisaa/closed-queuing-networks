% -------------------------------------------------------------------------
% Функция расчета
% -------------------------------------------------------------------------

function [lambda, j, t, to, V, f, No] = calculate(M, N, K, Q, MU, CS, errorRate)

    W = getW(); % матрица вероятностей нахождения заявки в текущем узле

    mostLoadedNode = getMostLoadedNode(); % наиболее загруженный узел

    l = 0;
    x = [];
    y = [];
    
    Bl = K(mostLoadedNode) * (MU(mostLoadedNode) / W(mostLoadedNode)) * (1 - 1 / N);
    
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

    % матрица вероятностей нахождения заявки в текущем узле
    
    function W = getW() 
        flag = false;
        
        for i = 1:M
            if (~ismember(1, Q(:,i)))
                flag = true;
                break;
            end
        end
        
        if flag
            A = transpose(Q);
            B = zeros(M,1);
        
            for i = 1:M
                A(i,i) = A(i,i) - 1;
            end
        
            A(M,:) = ones(1,M);
            B(M) = 1;
        
            W = linsolve(A, B);
        else
            W = ones(1,M);
        end
    end
        
    % Поиск самого загруженного узла
    
    function mostLoadedNode = getMostLoadedNode()
        max = 0;
        mostLoadedNode = 0;
        
        for i = 1:M
            value = W(i) / (K(i) * MU(i));
            if (value > max)
                max = value;
                mostLoadedNode = i;
            end
        end
    end
    
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
            Xi = zeros(M,1);
            F = zeros(M,1);
            alpha = zeros(M,1);
            beta = zeros(M,1);
            p_ij = zeros(M,M);
            
            for i = 1:M
                tmp = 0;
            
                for j = 1:M
                    tmp = tmp + (lambda_ij(j, i) / lambda(i))^2;
                    p_ij(j, i) = lambda_ij(j, i) / lambda(i);
                end
            
                Vi(i) = 1 / tmp;
                Wi(i) = 1 / (1 + 4 * (1 - ro(i))^2 * (Vi(i) - 1));
                Xi(i) = 1 + (max(CS(i), 0.2) - 1) / sqrt(K(i));
            end
            
            for i = 1:M
                sum_a = 0;
                for j = 1:M
                    beta(j,i) = Wi(i) * Q(j, i) * p_ij(j, i) * (1 - ro(j)^2);
                    sum_a = sum_a + p_ij(j, i) * ((1 - Q(j, i)) + Q(j, i) * ro(j)^2 * Xi(j));
                end
                alpha(i) = 1 + Wi(i) * (sum_a - 1);
            end
                
            for i = 1:M
                tmp = 0;
                for j = 1:M
                    tmp = tmp + ca(j) * beta(j, i);
                end
                 F(i) = ca(i) - alpha(i) - tmp;
            end
        
        end
        
        ca = fsolve(@caFun, ca0(M));
    end
    
    % P(0)

    function P0 = getP0(ro, m, factorial)
        P0 = 1 + ro^(m + 1) / (factorial(m) * (m - ro));
        
        for i = 1:m
            P0 = P0 + (ro^i / factorial(i));
        end
        
        P0 = 1 / P0;
    end
    
    % to без корректирующего фактора
    
    function toWCF = getToWCF(factorial, ro, ca)
        toWCF = zeros(M,1);
        
        for i = 1:M
            ro_i = ro(i);
            ca_i = ca(i);
            m_i = K(i);
            P0 = getP0(ro_i, m_i, factorial);
        
            to = (CS(i) + ca_i) / 2 * P0 * ro_i^m_i / (m_i * ((1 - ro_i / m_i)^2) * factorial(m_i) * MU(i));
            toWCF(i) = to;
        end
    end
    
    % Одна итерация алгоритма
    
    function [lambda, j, t, to] = doIteration(Bl)
        factorial = []; % массив факториалов
        kMax = max(K);
        
        fact = 1;
        for m = 1:(kMax + 1)
            fact = fact * m;
            factorial = [factorial fact];
        end

        lambda = Bl * W;
        lambda_ij = getLambdaij(lambda);
        [ca, ro] = getCai(lambda, lambda_ij);
        toWCF = getToWCF(factorial, ro, ca);
        
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

