% -------------------------------------------------------------------------
% Функция расчета
% -------------------------------------------------------------------------

function [lambda, j, t, to, V, f, No] = calculate(M, N, K, Q, MU, errorRate)

    W = getW(); % матрица вероятностей нахождения заявки в текущем узле

    mostLoadedNode = getMostLoadedNode(); % наиболее загруженный узел

    l = 0;
    x = [];
    y = [];
    
    Bl = (K(mostLoadedNode) / W(mostLoadedNode)) * MU(mostLoadedNode) * (1 - 1 / N);
    
    x = [x Bl];
    y = [y Bl];
    
    L = 0;
    
    while abs(L - N) >= errorRate
        l = l + 1;
        [lambda, j, No] = doIteration(Bl);
        j = transpose(j);
        L = sum(j);
        Bl = Bl * (N / L);
        y = [y Bl];
        if l > 1
            Bl = getBlWegstein(l, x, y);
        end
        x = [x Bl];
    end

    t = getT(j, lambda);
    to = getTw(t);
    V = getV(lambda);
    f = getF(V);

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
    
    % Одна итерация алгоритма
    
    function [lambda, j, No] = doIteration(Bl)
        lambda = Bl * W;

        j = [];
        No = [];

        F = [];
        kMax = max(K);
    
        factorial = 1;
        for m = 1:(kMax + 1)
            factorial = factorial * m;
            F = [F factorial];
        end
        
        for i = 1:M
            k = K(i);

            ro = lambda(i) / MU(i);
            sumH = 0;
        
            for s = 1:(k + 1)
                sumH = sumH + (ro^s) / F(s);
            end
        
            p0 = 1 / (1 + sumH + (ro^(k + 1)) / (F(k) * (k - ro)));
            No_i = (ro^(k + 1)) * p0 / (k * F(k) * ((1 - ro / k)^2));
            n_i = No_i + ro;
            j = [j n_i];
            No = [No No_i];
        end
    end
    
    % Метод Вегстейна
    
    function Bl_W = getBlWegstein(l, x, y)
        xk = x((l - 1):l);
        yk = y(l:(l + 1));
        
        Bl_W = yk(2) - ((yk(2) - yk(1)) * (yk(2) - xk(2))) / (yk(2) - yk(1) - xk(2) + xk(1));
    end
    
    % t

    function t = getT(j, lambda)
        t = zeros(1,M);
        
        for i = 1:M
            t(i) = j(i) / lambda(i);
        end
    end
    
    % t_ож
    
    function to = getTw(t)
        to = zeros(1,M);
        
        for i = 1:M
            to(i) = t(i) - 1 / MU(i);
        end
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

end

