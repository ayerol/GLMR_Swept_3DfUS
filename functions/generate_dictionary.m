function b = generate_dictionary(Fs,basis_shift,num_basis)

u = (0:1/Fs:8)';
construct_hrf = @(z) (z(2)*(((z(3)^z(1))*(u.^(z(1)-1)).*exp(-z(3)*u))/...
    gamma(z(1))));
b = [];
init_params = [3 1 4];

for idx_b = 1:num_basis
    hrf_param = shiftlat(init_params(1,:),(idx_b-1)*basis_shift/Fs);      
    b = cat(2,b,(construct_hrf(hrf_param)));
end

b(isnan(b)) = 0;


function [ phrf_new ] = shiftlat( phrf_orig , dt )

% SHIFTLAT transforms original HRF parameters into a new set of  parameters 
% that corresponds to a shift of the peak latencies with dt seconds.

phrf_new = phrf_orig;
phrf_new(1) = phrf_new(1) + dt*phrf_orig(3);
 
end

end