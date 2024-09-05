function [NS] = nash_sutcliffe(obs,mod)

% NaN values in the observations are replaced with empty cells
[r c] = find(isnan(obs));
obs(r) = [];
mod(r) = [];

numerator = sum((obs - mod).^2);
denom = sum((obs - mean(obs)).^2);

NS = 1 - (numerator/denom);

end

