% -------------------------------------------------------------------------
% Функция расчета
% -------------------------------------------------------------------------

function [lambda, j, t, to, V, f, No] = calculate(M, N, Q, R, MU, errorRate)

    W = getW(); % матрица вероятностей нахождения заявки в текущем узле

    mostLoadedNodes = getMostLoadedNodes(); % наиболее загруженный узел

    l = 0;
    Bl = zeros(1,R);
    
    for r = 1:R
        loadedNode = mostLoadedNodes(r);
        Bl(r) = (MU(r,loadedNode) / W(r, loadedNode)) * (1 - 1 / N(r));
    end
    
    x = Bl;
    y = Bl;
    
    error = errorRate * R;
    
    L = zeros(R,1);
    
    while sum(abs(L - N)) >= error
        l = l + 1;
        [to, lambda, j] = doIteration(Bl);
        j = transpose(j);
        L = sum(j, 2);
        for r = 1:R
            Bl(r) = Bl(r) * (N(r) / L(r));
        end
        y = [y; Bl];
        if l > 1
            Bl = getBlWegstein(l, x, y);
        end

        for r = 1:R
            if Bl(r) <= 0
                break;
            end
        end

        x = [x; Bl];
    end

    t = getT(j, lambda);
    V = getV(lambda);
    f = getF(V);
    No = getNo(j, to, t);

% -------------------------------------------------------------------------
% Вспомогательные функции
% -------------------------------------------------------------------------

    % матрица вероятностей нахождения заявки в текущем узле
    
    function W = getW() 
        W = zeros(R, M);

        for r = 1:R
            flag = false;
            Qr = Q(:,:,r);
        
            for i = 1:M
                if (~ismember(1, Qr(:,i)))
                    flag = true;
                    break;
                end
            end
        
            if flag
                A = transpose(Qr);
                B = zeros(M,1);
        
                for i = 1:M
                    A(i,i) = A(i,i) - 1;
                end
        
                A(M,:) = ones(1,M);
                B(M) = 1;
        
                W(r,:) = linsolve(A, B);
            else
                W(r,:) = ones(1,M);
            end
        end
    end
        
    % Поиск самых загруженных узлов

    function mostLoadedNodes = getMostLoadedNodes()
        mostLoadedNodes = zeros(1,R);
        
        for r = 1:R
            max = 0;
            maxIndex = 0;
        
            for i = 1:M
                value = W(r,i) / MU(r,i);
                if (value > max)
                    max = value;
                    maxIndex = i;
                end
            end
        
            mostLoadedNodes(r) = maxIndex;
        end
    end
    
    % Одна итерация алгоритма
    
    function [to, lambda, j] = doIteration(Bl)
        lambda = zeros(M,R);

        for i = 1:M 
            for r = 1:R
                lambda(i,r) = (Bl(r) * W(r,i));
            end
        end
        
        lambda_A = sum(lambda, 2);
        
        tau_A_r = zeros(M,R);
        for i = 1:M
            for r = 1:R
                tau_A_r(i,r) = (1 / MU(r,i)) * (lambda(i,r) / lambda_A(i));
            end
        end

        to = zeros(M,1);
        
        for i = 1:M
            tmp = zeros(1,R);
            for r = 1:R 
                tmp(r) = 1 / MU(r,i) * lambda(i,r) / lambda_A(i);
            end
            tmp_sum = sum(tmp);
            to(i) = lambda_A(i) * (tmp_sum^2) / (1 - lambda_A(i) * tmp_sum);
        end
        
        j = zeros(M,R);
        
        for i = 1:M
            ji = zeros(1,R);
            for r = 1:R
                ji(r) = lambda(i,r) * (to(i) + 1 / MU(r,i));
            end
            j(i,:) = ji;
        end
    end
    
    % Метод Вегстейна
    
    function Bl_W = getBlWegstein(l, x, y)
        xk = x((l - 1):l,:);
        yk = y(l:(l + 1),:);
        
        Bl_W = zeros(1,R);
        
        for r = 1:R
            Bl_W(r) = yk(2,r) - ((yk(2,r) - yk(1,r)) * (yk(2,r) - xk(2,r))) / (yk(2,r) - yk(1,r) - xk(2,r) + xk(1,r));
        end
    end

    % t

    function t = getT(j, lambda)
        t = zeros(R,M);
        
        for r = 1:R
            for i = 1:M
                t(r,i) = j(r,i) / lambda(i,r);
            end
        end
    end
    
    % Время цикла заявок разных классов для каждого из узлов сети

    function V = getV(lambda)
        V = zeros(R,M);
        
        for r = 1:R
            for i = 1:M
                V(r,i) = N(r) / lambda(i,r);
            end
        end
    end
    
    % Пропускная способность сети
    
    function f = getF(V)
        f = zeros(1,R);
        
        for r = 1:R
            f(r) = M / sum(V(r,:));
        end
    end
    
    % Число заявок в очереди
    
    function No = getNo(j, to, t)
        No = zeros(1,M);
        
        for i = 1:M
            Nq_i = zeros(1, R);
            for r = 1:R
                Nq_i(r) = j(r, i) * to(i) / t(r, i);
            end
            No(i) = sum(Nq_i);
        end
    end
end