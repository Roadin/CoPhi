function [best_norm, best_E, best_G] = find_best(best_norm, best_E, best_G,current_norm, current_E, current_G)
    for i=1:size(current_norm,2)
        if best_norm == 0 | best_norm > current_norm(i)
            best_norm = current_norm(i);
            best_E = current_E{i};
            best_G = current_G{i};
        end
    end
end
