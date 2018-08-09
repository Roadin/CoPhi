function [X_full_observe, F_full_observe, id_full_oberve, patient_full_observe]= delete_NP (F_org, X, id, patient)
    X_full_observe = zeros(size(X));
    F_full_observe = zeros(size(F_org));
    patient_full_observe = cell(size(patient));
    id_full_observe = zeros(size(id));

    rows = size(F_org, 1);
    counter = 1;
    for i=1:rows
        if F_org(i,1) ~= 0 && F_org(i,12) ~= 0
            X_full_observe(counter, :) = X(i,:);
            F_full_observe(counter, :) = F_org(i,:);
            id_full_observe(counter, :) = id(i,:);
            patient_full_observe(counter, :) = patient(i,:);
            counter = counter +1;
        end
    end
    X_full_observe(counter:end, :) = [];
    F_full_observe(counter:end, :) = [];
    id_full_observe(counter:end, :) = [];
    patient_full_observe(counter:end, :) = [];
end