function Tm = generate_measurement_tensor(T)

Tm=nan(size(T));
Ny = size(T,3);

round = 0;

for i=0:size(T,4)-1

    if round == 0 % for odd passovers, tensor is filled as in the normal
        % order of the slices

        Tm(:,:,mod(i,Ny)+1,i+1) = T(:,:,mod(i,Ny)+1,i+1);

        if i ~= 1 && mod(i,Ny) == Ny-1

            round = 1;

        end

    elseif round == 1 % for even passovers, tensor is filled in reverse
        % order of the slices (due to the continuous motion of the probe
        % while travelling back and forth)

        Tm(:,:,Ny-mod(i,Ny),i+1) = T(:,:,Ny-mod(i,Ny),i+1);

        if mod(i,Ny) == Ny-1

            round = 0;

        end

    end

end

end