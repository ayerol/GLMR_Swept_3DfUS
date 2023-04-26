%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Author: Aybuke Erol (a.erol@tudelft.nl)

% Date: 25.04.2023


% References

% [1] Michael Grant and Stephen Boyd. CVX: Matlab software for disciplined 
% convex programming, version 2.0 beta. http://cvxr.com/cvx, September 2013

% [2] B. Generowicz, S. Dijkhuizen, C. De Zeeuw, S. Koekkoek and P. 
% Kruizinga, "3D Functional Ultrasound using a Continuously Moving Linear 
% Stage," 2022 IEEE International Ultrasonics Symposium (IUS), 2022, 
% pp. 1-3, doi: 10.1109/IUS54386.2022.9957564.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all;
clear; clc;
rng default;

% Load the real imaging times and stimulus onsets from the fUS experiment
load('../data/tarr.mat'); % matrix of size Ny x T (# slices x time points)
load('../data/stim_on_off_times.mat'); % in terms of milliseconds

% Add the functions to search path
addpath(genpath('../functions/'))


%% %%%%%%%%%%%%%%%%%%%%%%%%%%% User Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Select a set of slices (so that the slice timings are realistic)
nys = 15:30:600; % to set a unified time-axis for reconstruction, use the
% central slice from each batch (assuming a batch size of 30)

Ny = length(nys); % number of slices

Nz = 20; % number of voxels along z-axis
Nx = 20; % number of voxels along x-axis

Nt = size(tarr,2); % number of time samples
% (i.e., number of times each slice was imaged)

rank = 2; % desired rank for decomposition
num_ROIs = 2; % number of regions of interest (ROIs)


%% %%%%%%%%%%%%%%%%%%%%%% Simulate ROI Activities %%%%%%%%%%%%%%%%%%%%%%%%%


% First create the time axis of the experiment by uniting the time samples
% from all considered slices

t_axis = [];
for nt = 1:Nt
    if mod(nt,2) == 1
        t_axis = cat(1,t_axis,tarr(nys,nt));
    else % the slice order is reversed during the returning of the probe
        t_axis = cat(1,t_axis,flipud(tarr(nys,nt)));
    end
end

% Create a time axis with uniform resolution to perform convolution for the
% design matrix

Fs = 100;
sim_t_axis = 0:1/Fs:180; % simulated time axis with uniform resolution

% Create the stimulus time course from the true stimulus onsets
% ( which will be used as input for the convolution )

stim = zeros(1,length(sim_t_axis));
for i = 1:2:length(stim_times)
    start = find(sim_t_axis == stim_times(i));
    endd = find(sim_t_axis == stim_times(i+1));
    stim(start:endd-1) = 1;
end
figure; plot(sim_t_axis,stim); xlim([0 sim_t_axis(end)]);
title('Stimulus');

% Perform convolution

y = cell(1,num_ROIs);
y_interp = cell(1,num_ROIs);

cmap = lines(7);
plot_size = 1/num_ROIs - 0.1;

for i = 1:num_ROIs

    [y{i},hrf] = generate_roi_response(stim,Fs);
    y_interp{i} = interp1(sim_t_axis,y{i},t_axis);
    if i == 1
        figure;
        plot(sim_t_axis,y{i},'LineWidth',2.5,'Color','k');
        hold on; plot(t_axis,y_interp{i},'*','Color',cmap(i,:));
        xlabel('Time (s)');
        set(gca,'FontSize',15); xlim([0 sim_t_axis(end)]);
        grid on; ylim([0 10]); yticks([0 4 8]);
        legend('True ROI Response','Sampled Points (from all slices)')
    end
end

figure; plot(0:1/Fs:8,hrf,'LineWidth',1.5);
xlabel('Time (s)'); title('Simulated HRF');


%% %%%%%%%%%%%%%%%%%%%%%%%% Simulate the Brain %%%%%%%%%%%%%%%%%%%%%%%%%%%%


noise = randn(length(t_axis),1);
noise(800:end) = noise(800:end) + 1.5; % add motion
figure; plot(noise); xlim([1 length(noise)]);
title('Common Background Activity');

% Create a brain using the common activity with varying amplitude per voxel

T = zeros([Nz Nx Ny length(t_axis)]);
br = rand(Nz,Nx,Ny); % brain map
for i = 1:length(noise)
    T(:,:,:,i) = br.*noise(i);
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%% Add the ROIs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


br_vis = zeros(size(br)); % only used for visualization
for i = 1:num_ROIs
    if i == 1
        roi_size = [5,5,9]; % depth, width, height
        roi_start = [2,2,1];
    elseif i == 2
        roi_size = [5,5,9]; % depth, width, height
        roi_start = [12,12,12];
    else
        roi_size = [randi([2,Nz-1]),randi([2,Nx-1]),randi([2,Ny])];
        roi_start = [randi([1,5]),randi([1,5]),randi([1,5])];
    end
    [br,T] = add_roi(br,roi_size,roi_start,T,y_interp{i});
    [br_vis,~] = add_roi(br_vis,roi_size,roi_start,T,y_interp{i});
end

figure;
h = slice(permute(br_vis,[1 3 2]),1:Ny,[],[]);
set(h,'EdgeColor','none',...
    'FaceColor','interp',...
    'FaceAlpha','interp')
alpha('color')
xlabel('Slice index')
zlabel('Voxel index in depth')
ylabel('Voxel index in width')
set(gca,'FontSize',16);
alphamap('increase',.05)
colormap parula;
cb = colorbar;
cbtitle = get(cb,'Title');
set(cbtitle,'String','Strength of the Response to Task');
set(cbtitle,'Rotation',270);
set(cbtitle,'FontSize',16.5);
set(cbtitle,'Position',[50 175 0]);
set(gca,'fontname','times');


%% %%%%%%%%%%%%%%% Create Data Array with Missing Values %%%%%%%%%%%%%%%%%%


Tm = generate_measurement_tensor(T);
Mm = reshape(Tm,Nz*Nx*Ny,[]); % matricise the data array


%% %%%%%%%%%%%%%%%%%%%%%%%%%% Visualization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Available samples at 2 arbitrary voxels

vx = [12,12,13];
f = figure;
p1 = plot(t_axis(find(~isnan(squeeze(Tm(vx(1),vx(2),vx(3),:))))),...
    squeeze(Tm(vx(1),vx(2),vx(3),...
    ~isnan(squeeze(Tm(vx(1),vx(2),vx(3),:))))),'o','Color',cmap(4,:));
p1.MarkerFaceColor = cmap(4,:);
p1.MarkerSize = 5;
hold on
vx = [3,18,1];
p2 = plot(t_axis(find(~isnan(squeeze(Tm(vx(1),vx(2),vx(3),:))))),...
    squeeze(Tm(vx(1),vx(2),vx(3),...
    ~isnan(squeeze(Tm(vx(1),vx(2),vx(3),:))))),'*','Color',cmap(5,:));
hold on
xlabel('Time (s)');
xticks(0:30:180)
yticks(-2:2:4)
ymin = -3; ymax = 5;
for i = 1:6
    p3 = fill([stim_times(i*2-1) stim_times(i*2-1) ...
        stim_times(i*2) stim_times(i*2)],...
        [ymin ymax ymax ymin],'r','HandleVisibility','off');
    set(p3,'facealpha',.1);
    set(p3,'EdgeColor','none');
end
ylim([ymin ymax])
legend([p1,p2,p3],{'In ROI';'Outside ROI';'Stimulus'},'Location',...
    'southeast','FontSize',14,'Orientation','horizontal');
set(gca,'FontSize',14);
set(gca,'fontname','times');
ax = gca;
f.Position = [440 377 560 200];


%% %%%%%%%%%%%%%%%%%%%%%%%%%% Regularizations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% General Linear Model (GLM)

basis_shift = 25; num_basis = 10; % shift the HRF parameters for every new
% column of the design matrix
bb = generate_dictionary(Fs,basis_shift,num_basis); % bb = hrf';
cmap = parula(num_basis); figure;
for i = 1:num_basis
    plot(0:1/Fs:8,bb(:,i)/max(bb(:,i)),'LineWidth',1.5,'Color',...
        cmap(num_basis-i+1,:));
    hold on;
end
hold off; set(gca,'YTick',[]); set(gca,'FontSize',12); xlabel('Time (s)');
title('Basis HRF Shapes');

design_matrix = zeros(length(t_axis),num_basis);
for i = 1:num_basis
    temp = conv(stim,bb(:,i));
    temp = temp(1:length(stim));
    design_matrix(:,i) = interp1(sim_t_axis,temp,t_axis,'previous');
end
design_matrix(isnan(design_matrix)) = 0;
figure; plot(t_axis,design_matrix,'LineWidth',1.5); set(gca,'FontSize',12); 
xlabel('Time (s)'); title('Design Matrix');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%% Parameters %%%%%
lambda_1 = .5; % GLM regularization
lambda_2 = 1e5; % sparsity of GLM coefficients
epsilon = 1e-3; % iteration-stopping criteria
maxit = 5; % max number of iterations if the algorithm doesnot converge

%%%%% Start the Algorithm %%%%%
idx = ~isnan(Mm); % observed indices
G = zeros(size(design_matrix,2),1);
G(randi([1,10])) = 1; % randomly initialize GLM coefficients
V = [design_matrix*G rand(size(Mm,2),1)]; % start with time

for it = 1:maxit

    if it > 1
        U_prev = U;
        V_prev = V;
    end

    cvx_begin quiet
    variable U(size(Mm,1),rank)
    minimize(myfunc(Mm,U,V,idx));
    cvx_end

    cvx_begin quiet
    variable V(size(Mm,2),rank)
    minimize(myfunc(Mm,U,V,idx) + lambda_1 * ...
        frosquared(V(:,1)-design_matrix*G));
    cvx_end

    cvx_begin quiet
    variable G(num_basis)
    minimize(lambda_1 * frosquared(V(:,1)-design_matrix*G) + ...
        lambda_2 * norm(G,1));
    cvx_end

    if it > 1

        delta_u = frosquared(U_prev-U)/frosquared(U);
        delta_v = frosquared(V_prev-V)/frosquared(V);

        if delta_u + delta_v < epsilon
            break;
        end

    end
    
    % % Uncomment below to visualize intermediate iterations
    % for r = 1:rank
    %     figure;
    %     h = slice(permute(reshape(U(:,r),Nz,Nx,Ny),[1 3 2]),1:Ny,[],[]);
    %     set(h,'EdgeColor','none',...
    %         'FaceColor','interp',...
    %         'FaceAlpha','interp')
    %     alpha('color')
    % 
    %     figure; plot(V(:,r));
    % end

end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


visualize_results(U,V,t_axis,stim_times,Nz,Nx,Ny);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function cost = myfunc(Mm,U,V,idx)
rec_Mm = U*V';
cost = sum(sum_square_abs(rec_Mm(idx)-Mm(idx)));
end

function cost = frosquared(Y)
cost = sum(sum_square_abs(Y));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%