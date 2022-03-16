
function y = SAcosPDF(data,kappa_m,capacity,kappa_r,kappa_f,s)

% Convert kappa into sd
J1=exp(kappa_m+kappa_f*cos(2*data.sample(:));

% First compute the number of items that get at least one slot
numRepresented = min(capacity, data.n(:));

% ... which we can use to compute the guess rate
g = 1 - numRepresented ./ data.n(:);

% Then pass around the slots evenly and compute the sd
slotsPerItemEvenly = floor(capacity ./ data.n(:));
worseSD = sd ./ max(sqrt(slotsPerItemEvenly),1); % to avoid Infs

% Count the items that get an extra slot and the resulting sd
numItemsWithExtraSlot = mod(capacity, data.n(:));
pExtraSlot = numItemsWithExtraSlot ./ numRepresented;
betterSD = sd ./ sqrt(slotsPerItemEvenly+1);

% Finally, compute probability
y =  (g(:)).*(1/360) + ...
    (1-g(:)).*(pExtraSlot(:)) .* vonmisespdf(data.errors(:),0,deg2k(betterSD(:))) + ...
    (1-g(:)).*(1-pExtraSlot(:)).*vonmisespdf(data.errors(:),0,deg2k(worseSD(:)));
end
