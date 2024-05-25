% матрица вероятностей нахождения заявки в текущем узле
    
function W = getW(M, Q) 
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
