%% simple FDTD solver for the wave equation

clear;
close all;
clc;

%% Choose solver
do_7_pt_stencil = 1;

id_sim = 'sim';

%% Physical/room parameters
L = 1;              % 1 meter
c = 0.5;            % speed of sound
duration = 0.75*L/c;    % in seconds
SR = 200;         % sample rate

eps = 1e-10;
lambda = 1/sqrt(3)-eps; % CFL condition ensured
k = 1/SR;               % time step
h = c*k/lambda;         % space step
N = floor(L/h);         % nb de pts spatiaux
h = L/N;
lambda = c*k/h;
nb_iteration = floor(duration/k);

%% Misc Study Parameters
plotOn = 1;
X = linspace(0,1,N);
save_gif = 1;
use_mult_scaling_factor = 0; % for plotting
mult_factor = 0.5;           % for plotting

%% Prepare data sets: first, specify the sampling frequency for sensors
f_sensors = 50;           % in Hz
DT = floor(SR/f_sensors); % every DT numerical iterations
idt = 1:DT:nb_iteration;
n_t = numel(idt);
f_sensors = 1/(DT*k);
n_sensors = 30;
rng default
n_exp = 1;
loc_sensors_full = zeros(n_sensors,3,n_exp);
for i=1:n_exp
    fprintf('Generating sensors %s out of %s \n',num2str(i),num2str(n_exp));
    loc_sensors = lhsdesign(n_sensors,3,'Criterion','maximin','Iterations',300);
    loc_sensors = 0.7*loc_sensors + 0.15;
    loc_sensors_full(:,:,i) = loc_sensors;
end

idspace = round(N*loc_sensors_full);

data_set = zeros(n_sensors*n_t,5,n_exp);
loc_sensors_toplot = loc_sensors(1:4,:);
figure(10)
plot3(loc_sensors_toplot(:,1),loc_sensors_toplot(:,2),loc_sensors_toplot(:,3),'x','color','r','MarkerSize',10,'linewidth',2);
title('Location of the sensors')
axis equal
xlim([0 1])
ylim([0 1])
zlim([0 1])

%% 1.2. Choose boundary conditions
is_neumann = 1;
is_EM_1 = 1 - is_neumann;

%% 1.3. Choose initial conditions
initial_cond.pos.is_excitation_initial = 1;
initial_cond.spd.is_excitation_initial = 0;
nr = 200;
uplot = zeros(nr,1);
vplot = zeros(nr,1);
rlist_u = [];
rlist_v = [];
x0_spd = [0 0 0];
R_spd = 0;
x0_pos = [0 0 0];
R_pos = 0;
u_surf = zeros(N,N);
v_surf = zeros(N,N);
%% 2.1. Generate Initial Conditions
u = zeros(N,N,N);

if initial_cond.pos.is_excitation_initial

    % for the cos SPEED ONLY :
%     initial_cond.x0 = [0.45 0.65 0.35];        % proportions of L
    % for the cos ring SPEED ONLY :
%     initial_cond.x0 = [0.3 0.6 0.7];
    % [0.25 0.35 0.5]
    % for the cos POSITION ONLY :
    % for the cos ring POSITION ONLY :


    initial_cond.pos.radius = 0.25;                % in proportions of L
%     dist_sensors = sqrt(sum((loc_sensors - initial_cond.pos.x0).^2,2));

    tell = 1;

    if tell == 1 %raised cosine
        % Initial condition parameters
%         initial_cond.amp = 10;
        % Amp = 100; for cos SPEED ONLY
        initial_cond.pos.x0 = [0.65 0.3 0.5];        % proportions of L
        Amp = 10; % WATCHOUT IS 5 FOR MIX;
        for i = 1:N
            for j = 1:N
                for z = 1:N
                    v = initial_cond.pos.x0 - [i/N j/N z/N];
                    if norm(v) < initial_cond.pos.radius
                        u(i,j,z) = 0.5*(1+cos( (pi/initial_cond.pos.radius)*norm(L.*v) ));
                    end
                end
            end
        end
        rmax = 2*initial_cond.pos.radius;
        rlist_u = linspace(0,rmax,nr);

        for i = 1:nr
            rh = rlist_u(i);
            if rh < initial_cond.pos.radius
                uplot(i) = Amp*0.5*(1+cos( (pi/initial_cond.pos.radius)*rh));
            end
        end

        u = Amp*u;
        integral = sum(sum(sum(u)))*h^3;
%         u = u/integral;
        tell_str = 'cos';
    elseif tell == 2 % radial gaussian centered at x_0
        sig = 0.01;
        for i = 1:N
            for j = 1:N
                for z = 1:N
                    v = initial_cond.pos.x0 - [i/N j/N z/N];
                    u(i,j,z) = exp(-0.5*norm(L*v)^2/sig^2);
                end
            end
        end
    elseif tell == 3 % centered indicator
        for i = 1:N
            for j = 1:N
                for z = 1:N
                    v = initial_cond.pos.x0 - [i/N j/N z/N];
                    if norm(v) < initial_cond.pos.radius
                        u(i,j,z) = 1;%*exp(-0.01/((i*k)^2));
                    end
                end
            end
        end
    elseif tell == 4 % Dirac

        id_x0 = round(L*N*initial_cond.pos.x0);
        u(id_x0(1),id_x0(2),id_x0(3)) = 1/h^3;

    elseif tell == 5 % Rectangle initial condition
        % in proportions of L
        xmin = 0.1;
        xmax = 0.9;
        ymin = 0.2;
        ymax = 0.3;
        zmin = 0.3;
        zmax = 0.4;
        x0 = L*[0.5*(xmin+xmax),0.5*(ymin+ymax),0.5*(zmin+zmax)];
        R = L*norm(x0-[xmin,ymin,zmin]);
        for i = 1:N
            xi = i/N;
            if (xmin < xi)&&(xi < xmax)
                for j = 1:N
                    yj = j/N;
                    if (ymin < yj)&&(yj < ymax)
                        for z = 1:N
                            zz = z/N;
                            if (ymin < yj)&&(yj < ymax)
                                u(i,j,z) = 1;%*exp(-0.01/((i*k)^2));
                            end
                        end
                    end
                end
            end
        end
    elseif tell == 6 % anneau cos

        initial_cond.pos.x0 = [0.5 0.5 0.5];        % proportions of L
        R1 = 0.15;
        R2 = 0.3;
        amp = 5;
        mid = 0.5*(R1+R2);
        diff = R2-R1;
        tell_str = 'cos_ring';
        for i = 1:N
            for j = 1:N
                for z = 1:N
                    v = initial_cond.pos.x0 - [i/N j/N z/N];
                    Nv = norm(v);
                    if (Nv < R2) && (Nv > R1)
                        u(i,j,z) = amp*(1+cos(2*pi*(Nv-mid)/diff));%*exp(-0.01/((i*k)^2));
                    end
                end
            end
        end
        rmax = 2*R2;
        rlist_u = linspace(0,rmax,nr);

        for i = 1:nr
            rh = rlist_u(i);
            if (rh < R2) && (rh > R1)
                uplot(i) = amp*(1+cos(2*pi*(rh-mid)/diff));
            end
        end

    end
    NZ = round(L*N*initial_cond.pos.x0(3));
    x0_pos = initial_cond.pos.x0;
    R_pos = initial_cond.pos.radius;
    u_surf = u(:,:,NZ);
    id_sim = [id_sim '_pos_' tell_str];
end


%%
v = zeros(N,N,N);

if initial_cond.spd.is_excitation_initial

    % for the cos SPEED ONLY :
    % for the cos ring SPEED ONLY :
    % [0.25 0.35 0.5]
    % for the cos POSITION ONLY :
%     initial_cond.spd.x0 = [0.65 0.3 0.5];        % proportions of L
    % for the cos ring POSITION ONLY :
%     initial_cond.spd.x0 = [0.5 0.5 0.5];        % proportions of L


    initial_cond.spd.radius = 0.15;                % in proportions of L
%     dist_sensors = sqrt(sum((loc_sensors - initial_cond.spd.x0).^2,2));

    tell = 1;

    if tell == 1 %raised cosine
        % Initial condition parameters
%         initial_cond.spd.amp = 10;
        % Amp = 100; for cos SPEED ONLY

        initial_cond.spd.x0 = [0.45 0.65 0.35];        % proportions of L
        Amp = 100; % or 200
        for i = 1:N
            for j = 1:N
                for z = 1:N
                    vh = initial_cond.spd.x0 - [i/N j/N z/N];
                    if norm(vh) < initial_cond.spd.radius
                        v(i,j,z) = 0.5*(1+cos( (pi/initial_cond.spd.radius)*norm(L.*vh) ));
                    end
                end
            end
        end
        v = Amp*v;
        integral = sum(sum(sum(v)))*h^3;
%         u = u/integral;
        rmax = 2*initial_cond.spd.radius;
        rlist_v = linspace(0,rmax,nr);

        for i = 1:nr
            rh = rlist_v(i);
            if rh < initial_cond.spd.radius
                vplot(i) = Amp*0.5*(1+cos( (pi/initial_cond.spd.radius)*rh));
            end
        end

    tell_str = 'cos';
    elseif tell == 2 % radial gaussian centered at x_0
        sig = 0.01;
        for i = 1:N
            for j = 1:N
                for z = 1:N
                    vh = initial_cond.spd.x0 - [i/N j/N z/N];
                    v(i,j,z) = exp(-0.5*norm(L*vh)^2/sig^2);
                end
            end
        end
    elseif tell == 3 % centered indicator
        for i = 1:N
            for j = 1:N
                for z = 1:N
                    vh = initial_cond.spd.x0 - [i/N j/N z/N];
                    if norm(vh) < initial_cond.spd.radius
                        v(i,j,z) = 1;%*exp(-0.01/((i*k)^2));
                    end
                end
            end
        end
    elseif tell == 4 % Dirac

        id_x0 = round(L*N*initial_cond.spd.x0);
        v(id_x0(1),id_x0(2),id_x0(3)) = 1/h^3;

    elseif tell == 5 % Rectangle initial condition
        % in proportions of L
        xmin = 0.1;
        xmax = 0.9;
        ymin = 0.2;
        ymax = 0.3;
        zmin = 0.3;
        zmax = 0.4;
        x0 = L*[0.5*(xmin+xmax),0.5*(ymin+ymax),0.5*(zmin+zmax)];
        R = L*norm(x0-[xmin,ymin,zmin]);
        for i = 1:N
            xi = i/N;
            if (xmin < xi)&&(xi < xmax)
                for j = 1:N
                    yj = j/N;
                    if (ymin < yj)&&(yj < ymax)
                        for z = 1:N
                            zz = z/N;
                            if (ymin < yj)&&(yj < ymax)
                                v(i,j,z) = 1;%*exp(-0.01/((i*k)^2));
                            end
                        end
                    end
                end
            end
        end
    elseif tell == 6 % anneau cos

        initial_cond.spd.x0 = [0.3 0.6 0.7];
        R1 = 0.05;
        R2 = 0.15; %R2 = 0.15;
        amp = 50;
        mid = 0.5*(R1+R2);
        diff = R2-R1;
        tell_str = 'cos_ring';
        for i = 1:N
            for j = 1:N
                for z = 1:N
                    vh = initial_cond.spd.x0 - [i/N j/N z/N];
                    Nvh = norm(vh);
                    if (Nvh < R2) && (Nvh > R1)
                        v(i,j,z) = amp*(1+cos(2*pi*(Nvh-mid)/diff));%*exp(-0.01/((i*k)^2));
                    end
                end
            end
        end
        rmax = 2*R2;
        rlist_v = linspace(0,rmax,nr);

        for i = 1:nr
            rh = rlist_v(i);
            if (rh < R2) && (rh > R1)
                vplot(i) = amp*(1+cos(2*pi*(rh-mid)/diff));
            end
        end
    end
    NZ = round(L*N*initial_cond.spd.x0(3));
    x0_spd = initial_cond.spd.x0;
    R_spd = initial_cond.spd.radius;
    v_surf = v(:,:,NZ);
    id_sim = [id_sim '_spd_' tell_str];
end

save('cond_ini.mat','u','v');

%% 2.3. Get Solver matrices

if do_7_pt_stencil
    if is_neumann
        Lpl_matrix_3D = kron_Lpl_3D_Neumann(N,N,N);
        B = (lambda^2).*Lpl_matrix_3D + 2*speye(N^3);
        A = -speye(N^3);
    elseif is_EM_1
        e = ones(N,1);
        Id_1D = speye(N);
        L1D = spdiags([e -2.*e e],-1:1,N,N);
        L1D(1,2) = 2;
        L1D(end,end-1) = 2;

        edges = sparse(N,N);
        edges(1,1) = 1;
        edges(end,end) = 1;
        B1D = (lambda^2).*L1D + 2*speye(N);
        Edges_3D = kron(edges,kron(Id_1D,Id_1D)) + kron(Id_1D,kron(edges,Id_1D)) + kron(Id_1D,kron(Id_1D,edges));
        C = speye(N^3) + lambda*Edges_3D;
        A = lambda*Edges_3D-speye(N^3);
        B = (lambda^2)*(kron(L1D,kron(Id_1D,Id_1D)) + kron(Id_1D,kron(L1D,Id_1D)) + kron(Id_1D,kron(Id_1D,L1D)) )  + 2*speye(N^3);

        B = C^(-1)*B;
        A = C^(-1)*A;
    end
end

%% 2.4. Launch solver

% Initialize study variables
u3 = zeros(N^3,1);
y = zeros(nb_iteration,1);
H = zeros(nb_iteration,1);       % total energy
T = zeros(nb_iteration,1);       % kinetic
V = zeros(nb_iteration,1);       % potential
T_1 = zeros(nb_iteration,1);     % other kinetic

% Reshape initial conditions to vector format
u1 = reshape(u,[N^3,1]);
u2 = reshape(u,[N^3,1]) + k*reshape(v,[N^3,1]);

Nsq = N*N;
Nm = N - 1;

lba = zeros(N^2,1);
lf = zeros(N^2,1);
lr = zeros(N^2,1);
ll = zeros(N^2,1);
lbo = zeros(N^2,1);
lt = zeros(N^2,1);
CfaceC = 1/(1 + lambda);
CfaceA = lambda - 1;
order2 = 1;

cnt_t = 1;
cnt_im = 1;
tic;
fprintf('Launching simulation... \n');

fig = figure(1);

for t = 1:nb_iteration

    % interior domain and boundary conditions
    u3 = B*u2 + A*u1;

    if ismember(t,idt)
        tmat = duration*t/nb_iteration;
        for id_exp =1:n_exp

            for id_sensor = 1:n_sensors
                loc_s_h = idspace(id_sensor,:,id_exp);
                id_full = loc_s_h(1) + N*(loc_s_h(2)-1) + N*N*(loc_s_h(3)-1);
                id_data = id_sensor + n_sensors*(cnt_t-1);
                data_set(id_data,:,id_exp) = [loc_sensors_full(id_sensor,:,id_exp),tmat,u3(id_full)];
            end
        end
        cnt_t = cnt_t+1;
    end

    if plotOn
        if rem(t,3) == 0
        uplot = reshape(u3,[N,N,N]);

        if 0
            %subplot(1,2,1)
    %         surf(log(abs(uplot(:,:, floor(N*sources.trajectory.curve(3)) ))) )
            surf(uplot(:,:, floor(N*initial_cond.x0(3)) ) )
            hold on
            plot3(sources.trajectory.curve2plot(2,:),sources.trajectory.curve2plot(1,:),sources.trajectory.curve2plot(3,:),'w');
            plot3(sources.trajectory.curve2plot(2,t),sources.trajectory.curve2plot(1,t),sources.trajectory.curve2plot(3,t),'Color','r','Marker','x','MarkerSize',15,'LineWidth',2);
    %         plot3([0,sources.trajectory.curve2plot(1,t)],[0,sources.trajectory.curve2plot(2,t)],[0,sources.trajectory.curve2plot(3,t)],'w','x');
            colorbar
            if use_mult_scaling_factor
                caxis(mult_factor.*[-1 1])
            end
            axis equal
            axis ij
            hold off
            shading interp
%             shading flat
            view(0,90)
            drawnow

            %subplot(1,2,2)
            %plot(vecOptiCoeffs(1:end-1)) %,'o','LineWidth',3
            %ylim([-0.5 1])
            %drawnow
        else
%             surf(log(abs(uplot(:,:, floor(N*sources.trajectory.curve(3)) ))) )
%      surf(uplot(:,:, floor(N*initial_cond.x0(3)))-uplot(:,:, floor(N*initial_cond.x0(3)))' )

            surf(X,X,uplot(:,:, NZ))
%             surf(reshape(uplot(:,floor(N*sources.trajectory.curve(2)),:),N,N) )
            hold on
            plot3(loc_sensors(:,1),loc_sensors(:,2),ones(n_sensors,1),'x','color','r','MarkerSize',10,'linewidth',2);

            colorbar
            if use_mult_scaling_factor
                caxis(mult_factor.*[-0.1 1])
            end
            axis square
            axis xy
            hold off
            shading interp
%             shading flat
            view(0,90)
            title(sprintf('t = %s out of %s',num2str(t*duration/nb_iteration),num2str(duration)))
            drawnow

            if save_gif
                frame = getframe(fig);
                im{cnt_im} = frame2im(frame);
                cnt_im = cnt_im+1;
            end
        end
        end
    end

    if do_energy_test

       % Reshape solution to matrix format
       u3mat = reshape(u3,N,N,N);
       u2mat = reshape(u2,N,N,N);

       % Compute
       xtdiff = u3mat(2:end,:,:)-u3mat(1:end-1,:,:);
       ytdiff = u3mat(:,2:end,:)-u3mat(:,1:end-1,:);
       ztdiff = u3mat(:,:,2:end)-u3mat(:,:,1:end-1);
       xtmdiff = u2mat(2:end,:,:)-u2mat(1:end-1,:,:);
       ytmdiff = u2mat(:,2:end,:)-u2mat(:,1:end-1,:);
       ztmdiff = u2mat(:,:,2:end)-u2mat(:,:,1:end-1);
       dotpx = xtdiff.*xtmdiff;
       dotpy = ytdiff.*ytmdiff;
       dotpz = ztdiff.*ztmdiff;

       T(t) = 0.5*h^3*sumsqr((1/k).*(u3-u2));                           % kinetic energy
       V(t) = 0.5*c^2*h*(sum(dotpx(:))+sum(dotpy(:))+sum(dotpz(:)));    % potential energy
       H(t) = T(t)+V(t);
    end

    % Update pressure fields
    u1 = u2;
    u2 = u3;

    % Show duration already simulated
    if rem(t,5) == 0
        timer = t*k;
        fprintf('%s s out of %s s || %s %% finished \n',num2str(timer),num2str(duration),num2str(round(1000*timer/duration)/10));
    end
end
t_matrix = toc;
fprintf('Simulation finished! \n');

%% save stuff
if save_gif
    fprintf('Saving gif...')
    tmat = now;
    tmat = floor(1e7*tmat);
    tmat = rem(tmat,1e5);
    duration_gif = 3; % 3 seconds
    filename = sprintf('%s_compact_wave_propag_MATLAB.gif',num2str(tmat));
    full_fname = [pwd sprintf('/gifs_images/%s',filename)];
    delay_t = duration_gif/(cnt_im-1);
    for idt_im = 1:(cnt_im-1)
        [A,map] = rgb2ind(im{idt_im},256);
        if idt_im == 1
            imwrite(A,map,full_fname,'gif','LoopCount',Inf,'DelayTime',delay_t);
        else
            imwrite(A,map,full_fname,'gif','WriteMode','append','DelayTime',delay_t);
        end
    end
    fprintf('Done! \n')
end

save('t_forsave.mat','n_t','n_sensors','c','mult_factor','loc_sensors') % ,'tmat'
save('physical_params.mat', 'x0_pos','R_pos','x0_spd','R_spd','c')
date = datetime('now');

fprintf('Saving data set output...');
reshape_data = 1;
if reshape_data % put in sort of time series format
    new_data = zeros(size(data_set));
    n_data = n_sensors*n_t;
    for id_exp=1:n_exp
        for i=1:n_sensors
            id_sensor_h = i:n_sensors:n_data;
            id_new_h = (1+(i-1)*n_t):(i*n_t);
            new_data(id_new_h,:,id_exp) = data_set(id_sensor_h,:,id_exp);
        end
    end

    noisy_full = new_data;

    for id_exp=1:n_exp
        writematrix(new_data(:,:,id_exp),sprintf('exp_%s_data_set_noiseless.csv',num2str(id_exp)));
        noise_lvl = 9e-2; % choose the suitable noise level
        noise = noise_lvl*randn(n_data,1);
        noisy_obs = new_data(:,:,id_exp);
        noisy_obs(:,end) = noisy_obs(:,end) + noise;
        noisy_full(:,5,id_exp) = noisy_obs(:,end);
        writematrix(noisy_obs,sprintf('exp_%s_data_set_noisy.csv',num2str(id_exp)));
    end
 %%
    xt = k*(idt-1);
    wc = 2*pi/(xt(2)-xt(1));
    omega = linspace(-wc/2,wc/2,n_t);
    for i=1:n_sensors
        id_new_h = (1+(i-1)*n_t):(i*n_t);
        s_h = new_data(id_new_h,5);
        fft_h = fftshift(fft(s_h));
        figure(10)
        subplot(3,1,1)
        plot(xt,s_h,'x')
        subplot(3,1,2)
        plot(omega,abs(fft_h))
        subplot(3,1,3)
        plot(omega,(180/pi)*unwrap(angle(fft_h)))

    end
    %% show signals
    for id_exp=1:n_exp
        figure(11)
        subplot(2,1,1)
        plot(new_data(:,5,id_exp))
        hold on
        for i=1:n_sensors
            xline(i*n_t)
        end
        hold off
        xlim([1 n_t*n_sensors])
        subplot(2,1,2)
        plot(noisy_full(:,5,id_exp))
        hold on
        for i=1:n_sensors
            xline(i*n_t)
        end
        hold off
        xlim([1 n_t*n_sensors])
    end
else
    writematrix(data_set,'data_set.csv');
end
fprintf('Done! \n')

data_show = [new_data(:,5,id_exp),noisy_full(:,5,id_exp)];
data_show_str = [id_sim '.mat'];
save(data_show_str,'data_show','rlist_u','rlist_v','uplot','vplot','u_surf','v_surf'); %,'u','v'
