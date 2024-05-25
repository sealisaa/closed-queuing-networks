% -------------------------------------------------------------------------
% Функция расчета
% -------------------------------------------------------------------------

function [lambda, j, t, to, V, f, No, MU_A, lambda_A, q, CS_A, ca, ro] = calculate(M, N, Q, R, MU, CS, errorRate)

    W = getW(); % матрица вероятностей нахождения заявки в текущем узле

    mostLoadedNodes = getMostLoadedNodes(); % наиболее загруженный узел

    l = 0;
    x = [];
    y = [];
    flag = 1;
    
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
        [lambda, j, t, to, MU_A, lambda_A, q, CS_A] = doIteration(Bl);
        j = transpose(j);
        L = sum(j, 2);
        for r = 1:R
            Bl(r) = Bl(r) * (N(r) / L(r));
        end
        y = [y; Bl];
        if l > 1
            [Bl, flag] = getBlWegstein(l, x, y);
        end
        x = [x; Bl];

        for r = 1:R
            if Bl(r) <= 0 
                flag = 1;
            end
        end

        if (flag == 0 || ismember(0, Bl))
            break;
        end
    end

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
    
    % lambda_ij
    
    function lambda_ij = getLambdaij(lambda)
       lambda_ij = zeros(M, M, R);
    
        for r = 1:R
            lambda_ij_r = zeros(M);
            for i = 1:M
                for j = 1:M
                    lambda_ij_r(i, j) = lambda(i, r) * Q(i, j, r);
                end
            end
            lambda_ij(:,:,r) = lambda_ij_r;
        end
    end
    
    % q
    
    function q = getQ(lambda_ij, lambda_A)
        for i = 1:M
            for j = 1:M
                tmp_sum = 0;

                for r = 1:R
                    tmp_sum = tmp_sum + lambda_ij(i,j,r);
                end

                q(i,j) = tmp_sum / lambda_A(i);
            end
        end
    end
    
    % ca0
    
    function [ca] = ca0(M)
        ca = zeros(M,1);
        ca(1)=1;
    end

    % ca

    function [ca, ro, CS_A] = getCai(q, lambda, lambda_ij, lambda_A, MU_A)
        lambda_Q = zeros(M);

        for i = 1:M
            for j = 1:M
                lambda_Q(i, j) = lambda_A(i) * q(i, j);
            end
        end
        
        CS_A = zeros(M,1);
        
        for i = 1:M
            tmp = 0;

            for r = 1:R
                tmp = tmp + (lambda(i, r) / lambda_A(i));
            end

            CS_A(i) = tmp * (CS(r, i) + 1) * ((MU_A(i) / MU(r, i))^2) - 1;
        end
    
        ro = zeros(M,1);
    
        for i = 1:M
            ro(i) = lambda_A(i) / MU_A(i);
        end
        
        function [F] = caFun(ca)
            Wi = zeros(M,1);
            Vi = zeros(M,1);
            F = zeros(M,1);
            
            for i = 1:M
                tmp = 0;
            
                for j = 1:M
                    tmp = tmp + (lambda_ij(j, i) / lambda_A(i))^2;
                end
            
                Vi(i) = 1 / tmp;
                Wi(i) = 1 / (1 + 4 * (1 - ro(i))^2 * (Vi(i) - 1));
            end
            
            for i = 1:M
                tmp = 0;
            
                for j = 1:M
                    ca_ij = q(j, i) * (ro(j)^2 * CS_A(j) + (1 - ro(j)^2) * ca(j)) + 1 - q(j,i);
                    tmp = tmp + lambda_ij(j, i) * ca_ij / lambda_A(i);
                end
                F(i) = Wi(i) * tmp + 1 - Wi(i) - ca(i);
            end
        
        end
        
        ca = fsolve(@caFun, ca0(M));
    end
    
    % to без корректирующего фактора
    
    function toWCF = getToWCF(ro, ca, MU_A, CS_A)
        toWCF = zeros(M,1);
        
        for i = 1:M
            ro_i = ro(i);
            ca_i = ca(i);
            to = ro_i * (CS_A(i) + ca_i) / (2 * MU_A(i) * (1 - ro_i));
            if ca_i < 1
                to = to * exp((ro_i - 1) * ((1 - ca_i)^2) / (1.5 * ro_i * (CS_A(i) + ca_i)));
            else
                to = to * exp((ro_i - 1) * (ca_i - 1) / (4 * CS_A(i) + ca_i));
            end
            toWCF(i) = to;
        end
    end
    
    % Одна итерация алгоритма
    
    function [lambda, j, t, to, MU_A, lambda_A, q, CS_A] = doIteration(Bl)
        lambda = zeros(M,R);
    
        for i = 1:M 
            for r = 1:R
                lambda(i,r) = (Bl(r) * W(r,i));
            end
        end
        
        lambda_A = sum(lambda, 2);
        lambda_ij = getLambdaij(lambda);
        
        q = getQ(lambda_ij, lambda_A);
    
        tau_A_r = zeros(M,R);
        for i = 1:M
            for r = 1:R
                tau_A_r(i,r) = (1 / MU(r,i)) * (lambda(i,r) / lambda_A(i));
            end
        end
        
        tau_A = sum(tau_A_r, 2);
        
        MU_A = zeros(M,1);
        
        for i = 1:M
            MU_A(i) = 1 / tau_A(i);
        end
    
        [ca, ro, CS_A] = getCai(q, lambda, lambda_ij, lambda_A, MU_A);
        toWCF = getToWCF(ro, ca, MU_A, CS_A);
    
        N_sum = sum(N);
        
        to = zeros(M,1);
        t = zeros(M,R);
        j = zeros(M,R);
        
        for i = 1:M
            to_i = toWCF(i) * (N_sum - 1) / (N_sum + toWCF(i) * MU_A(i));
            to(i) = to_i;
    
            for r = 1:R
                t(i,r) = to_i + 1 / MU(r, i);
                j(i, r) = t(i,r) * lambda(i,r);
            end
        end
    end
    
    % Метод Вегстейна
    
    function [Bl_W, flag] = getBlWegstein(l, x, y)
        flag = 1;
        xk = x((l - 1):l,:);
        yk = y(l:(l + 1),:);
        
        Bl_W = zeros(1,R);
        
        for r = 1:R
            tmp = (yk(2,r) - yk(1,r) - xk(2,r) + xk(1,r));

            if tmp == 0
                flag = 0;
                return;
            end

            Bl_W(r) = yk(2,r) - ((yk(2,r) - yk(1,r)) * (yk(2,r) - xk(2,r))) / tmp;
        end
    end
    
    % Время цикла заявок для каждого из узлов сети
    
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
                Nq_i(r) = j(r, i) * to(i) / t(i, r);
            end
            No(i) = sum(Nq_i);
        end
    end
end