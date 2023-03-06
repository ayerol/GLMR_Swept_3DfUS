function [y,hrf] = generate_roi_response(stim,Fs)

u = 0:1/Fs:8;
construct_hrf = @(z) (z(2)*(((z(3)^z(1))*(u.^(z(1)-1)).*exp(-z(3)*u))/...
    gamma(z(1)))); % gamma-model for generating the HRFs
hrf = construct_hrf([3 rand 4]); % values used in the paper, can be changed
y = conv(hrf,stim);
y = y(1:length(stim))/10;

end