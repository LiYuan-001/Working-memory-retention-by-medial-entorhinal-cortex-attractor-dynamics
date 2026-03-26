function remainInd = removeConsBin(vec,consBinThres)

% Count occurrences
[uniqueVals, ~, idx] = unique(vec);
counts = accumarray(idx, 1);

% Find values that appear more than consBinThres times
valsToRemove = uniqueVals(counts > consBinThres);

% Remove them from the vector
remainInd = ~ismember(vec, valsToRemove);

end