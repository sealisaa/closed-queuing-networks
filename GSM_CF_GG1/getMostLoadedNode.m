% Поиск самого загруженного узла
    
function mostLoadedNode = getMostLoadedNode(M, MU, W)
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