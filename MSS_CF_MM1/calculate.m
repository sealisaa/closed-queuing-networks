% -------------------------------------------------------------------------
% Функция расчета
% -------------------------------------------------------------------------

function [lambda, j, t, to, V, f, No] = calculate(M, N, Q, MU, errorRate)

    W = getW(); % матрица вероятностей нахождения заявки в текущем узле

    mostLoadedNode = getMostLoadedNode(); % наиболее загруженный узел

    l = 0;
    x = [];
    y = [];
    
    Bl = (MU(mostLoadedNode) / W(mostLoadedNode)) * (1 - 1 / N);
    
    x = [x Bl];
    y = [y Bl];
    
    L = 0;
    
    while abs(L - N) >= errorRate
        l = l + 1;
        [lambda, j] = doIteration(Bl);
        j = transpose(j);
        L = sum(j);
        Bl = Bl * (N / L);
        y = [y Bl];
        if l > 1
            Bl = getBlWegstein(l, x, y);
        end

        if Bl <= 0
            break;
        end

        x = [x Bl];
    end

    t = getT(j, lambda);
    to = getTo(t);
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
            value = W(i) / MU(i);
            if (value > max)
                max = value;
                mostLoadedNode = i;
            end
        end
    end
    
    % Одна итерация алгоритма
    
    function [lambda, j] = doIteration(Bl)
        lambda = Bl * W;

        j = zeros(M,1);
        
        for i = 1:M
            tmp = MU(i) * (1 - lambda(i) / MU(i));
            j(i) = lambda(i) * (((N - 1) / N) * (N / (N + lambda(i) / tmp)) * (lambda(i) / MU(i) / tmp) + 1 / MU(i));
        end
    end
    
    % Метод Вегстейна
    
    function Bl_W = getBlWegstein(l, x, y)
        xk = x((l - 1):l);
        yk = y(l:(l + 1));
        
        Bl_W = yk(2) - ((yk(2) - yk(1)) * (yk(2) - xk(2))) / ((yk(2) - yk(1)) - (xk(2) + xk(1)));
    end

    % t
    
    function t = getT(j, lambda)
        t = zeros(1,M);
        
        for i = 1:M
            t(i) = j(i) / lambda(i);
        end
    end
    
    % t_ож
    
    function to = getTo(t)
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
    
    % Число заявок в очереди
    
    function No = getNo(j, to, t)
        No = zeros(1,M);
    
        for i = 1:M
            No(i) = j(i) * to(i) / t(i);
        end
    end
end

