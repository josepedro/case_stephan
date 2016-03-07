
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This is a simple implementation of the LBGK model.
%% By Andrey R. da Silva, August 2010
%%
%% The code does not take into account eny specific boundaru condiition.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc
close all

tic

% setting size of pipe
radius_pipe = 160;
length_pipe = 160*14;
length_anechoic_condition = 30;

% Block 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Lattice size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nr = radius_pipe + 11;                    % Number of lines   (cells in the y direction)
Mc = length_pipe + 2*length_anechoic_condition + 2;                    % Number of columns (cells in the x direction)

% Block 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Physical parameters (macro)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c_p = 340;                    % Sound velocity on the fluid [m/s]
rho_p = 1.2;                  % physical density [kg/m^3]
rho_p = 1;                    % Fluid density  [kg/m^3]
Lx = .5;                      % Maximun dimenssion in the x direction [m]
Ly = 0.075;                  % Maximun dimenssion on th y direction  [m]
Dx = Ly/radius_pipe                    % Lattice space (pitch)
Dt = (1/sqrt(3))*Dx/c_p       % lattice time step

% Block 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Lattice parameters (micro - lattice unities)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
omega = 1.98;                                     % Relaxation frequency
tau = 1/omega;                                    % Relaxation time
rho_l = 1;                                        % avereged fluid density (latice density
cs = 1/sqrt(3);                                   % lattice speed of sound
cs2 = cs^2;                                       % Squared speed of sound cl^2
visc = cs2*(1/omega-0.5);                         % lattice viscosity
visc_phy = visc*(Dx^2)/Dt;                        % physical kinematic viscosity

% Block 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lattice properties for the D2Q9 model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_c = 9;                                             % number of directions of the D2Q9 model
C_x = [1 0 -1  0 1 -1 -1  1 0];                       % velocity vectors in x
C_y = [0 1  0 -1 1  1 -1 -1 0];                       % velocity vectors in y
w0=16/36. ; w1=4/36. ; w2=1/36.;                    % lattice weights
W = [w1 w1 w1 w1 w2 w2 w2 w2 w0];
f1=3.;
f2=4.5;
f3=1.5;                                             % coef. of the f equil.

% Array of distribution and relaxation functions
f=zeros(Nr,Mc,N_c);                                 
feq=zeros(Nr,Mc,N_c);

% Filling the initial distribution function (at t=0) with initial values
f(:,:,:)=rho_l/9;
ux = zeros(Nr, Mc);
uy = zeros(Nr, Mc);

% Calculando a condicao anecoica
% 4.0.1 - Adding conditions anechoic
% condicao anecoica para cima
growth_delta = 0.5;
size_up_anechoic = 10;
[sigma_mat9_cima Ft_cima] = build_anechoic_condition(Mc, ... 
Nr, size_up_anechoic, growth_delta);
% condicao anecoica para esquerda
growth_delta = -1;
[sigma_mat9_esquerda Ft_esquerda] = build_anechoic_condition(Mc, ... 
Nr, length_anechoic_condition, growth_delta);
% condicao anecoica para direita
growth_delta = 1;
[sigma_mat9_direito Ft_direito] = build_anechoic_condition(Mc, ... 
Nr, length_anechoic_condition, growth_delta);

% Vendo pontos de barreira da xirizuda
xl = [2 2 2300+1 2300+1];
yl = [1 160 160 1];
[vec1,vec2,vec3,vec4,vec5,vec6,vec7,vec8] = crossing3(Nr,Mc,xl,yl);
%xl = [250 250];
%yl = [1 50];
xl = [1 2300];%Lt=230 voxel; Ltubo=200 voxel
yl = [161 161];
[vec1_duto,vec2_duto,vec3_duto,vec4_duto,vec5_duto,vec6_duto,vec7_duto,vec8_duto] = ... 
crossing3(Nr,Mc,xl,yl);

%Block 522
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Begin the iteractive process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construindo chirp
a=yl(2)-1;
total_time = 25200; % meia hora = 20*Mc*sqrt(3)
final_time_chirp = total_time - 2*Mc*sqrt(3);
times = 0 : final_time_chirp;
initial_frequency = 0;
frequency_max_lattice = 4*cs/(2*pi*a);
source_chirp = chirp(times, ... 
initial_frequency, times(end), frequency_max_lattice);
% vetores de pressao e velocidade de particula
pressure1(1:total_time) = 0;
pressure2(1:total_time) = 0;
pressure3(1:total_time) = 0;
particle_velocity1(1:total_time) = 0;
particle_velocity2(1:total_time) = 0;
particle_velocity3(1:total_time) = 0;
for ta = 1 : total_time
    
    % Block 5.1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% propagation (streaming)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    f(:,:,1) = [f(:,1:2,1) f(:,2:Mc-1,1)];
    f(:,:,2) = [f(1:2,:,2);f(2:Nr-1,:,2)];
    f(:,:,3) = [f(:,2:Mc-1,3) f(:,Mc-1:Mc,3)];
    f(:,:,4) = [f(2:Nr-1,:,4);f(Nr-1:Nr,:,4)];
    f(:,:,5) = [f(:,1:2,5) f(:,2:Mc-1,5)];
    f(:,:,5) = [f(1:2,:,5);f(2:Nr-1,:,5)];
    f(:,:,6) = [f(:,2:Mc-1,6) f(:,Mc-1:Mc,6)];
    f(:,:,6) = [f(1:2,:,6);f(2:Nr-1,:,6)];
    f(:,:,7) = [f(:,2:Mc-1,7) f(:,Mc-1:Mc,7)];
    f(:,:,7) = [f(2:Nr-1,:,7);f(Nr-1:Nr,:,7)];
    f(:,:,8) = [f(:,1:2,8) f(:,2:Mc-1,8)];
    f(:,:,8) = [f(2:Nr-1,:,8);f(Nr-1:Nr,:,8)];

    G=f;
    f(vec1)=G(vec3);
    f(vec3)=G(vec1);
    f(vec2)=G(vec4);
    f(vec4)=G(vec2);
    f(vec5)=G(vec7);
    f(vec7)=G(vec5);
    f(vec6)=G(vec8);
    f(vec8)=G(vec6);
    
    G=f;
    f(vec1_duto)=G(vec3_duto);
    f(vec3_duto)=G(vec1_duto);
    f(vec2_duto)=G(vec4_duto);
    f(vec4_duto)=G(vec2_duto);
    f(vec5_duto)=G(vec7_duto);
    f(vec7_duto)=G(vec5_duto);
    f(vec6_duto)=G(vec8_duto);
    f(vec8_duto)=G(vec6_duto);

    % Block 5.2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % recalculating rho and u
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rho=sum(f,3); 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calculando uma fonte ABC dentro do duto
    % direita = 0.5
    if ta <= length(source_chirp)
        density_source = rho_l + 0.0001*source_chirp(ta);
        Ux_t = (0.0001*source_chirp(ta))/sqrt(3);
    else
        density_source = rho_l;
        Ux_t = 0;
    end
    Uy_t = 0;
    point_y = 1;
    distance_y = radius_pipe;
    point_x = 2;
    distance_x = 30;
    direction = 0.5;
    [sigma_source Ft_source] = build_source_anechoic(Nr, Mc, ...
    density_source, Ux_t, Uy_t, point_y, ...
    point_x, distance_x, distance_y, direction);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ta == 1
    %rho(15, 100) = rho_l + 0.0001 *sin(ta);
    end
    
    rt0= w0*rho;
    rt1= w1*rho;
    rt2= w2*rho;

    % Determining the velocities according to Eq.() (see slides)    
    ux = (C_x(1).*f(:,:,1)+C_x(2).*f(:,:,2)+C_x(3).*f(:,:,3)+C_x(4).*f(:,:,4)+C_x(5).*f(:,:,5)+C_x(6).*f(:,:,6)+C_x(7).*f(:,:,7)+C_x(8).*f(:,:,8))./rho ;
    uy = (C_y(1).*f(:,:,1)+C_y(2).*f(:,:,2)+C_y(3).*f(:,:,3)+C_y(4).*f(:,:,4)+C_y(5).*f(:,:,5)+C_y(6).*f(:,:,6)+C_y(7).*f(:,:,7)+C_y(8).*f(:,:,8))./rho ;
        

    % Block 5.3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Determining the relaxation functions for each direction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
    uxsq=ux.^2; 
    uysq=uy.^2; 
    usq=uxsq+uysq; 
    u=sqrt(usq);
    feq(:,:,1)= rt1 .*(1 +f1*ux +f2.*uxsq -f3*usq);
    feq(:,:,2)= rt1 .*(1 +f1*uy +f2*uysq -f3*usq);
    feq(:,:,3)= rt1 .*(1 -f1*ux +f2*uxsq -f3*usq);
    feq(:,:,4)= rt1 .*(1 -f1*uy +f2*uysq -f3*usq);
    feq(:,:,5)= rt2 .*(1 +f1*(+ux+uy) +f2*(+ux+uy).^2 -f3.*usq);
    feq(:,:,6)= rt2 .*(1 +f1*(-ux+uy) +f2*(-ux+uy).^2 -f3.*usq);
    feq(:,:,7)= rt2 .*(1 +f1*(-ux-uy) +f2*(-ux-uy).^2 -f3.*usq);
    feq(:,:,8)= rt2 .*(1 +f1*(+ux-uy) +f2*(+ux-uy).^2 -f3.*usq);
    feq(:,:,9)= rt0 .*(1 - f3*usq);
    
    % Block 5.4
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Collision (relaxation) step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Condicao Axissimetrica
        % termo de primeira ordem
            % y por x
    h_1_leaf = zeros(Nr, Mc);
    radius = 1:Nr;
    radius = radius';
    for column = 1:Mc
        h_1_leaf(:, column) = uy(:, column)./radius;
    end
            % construindo a matriz com os pesos de cada direcao
    h_1 = zeros(Nr,Mc,N_c);
    for direction = 1:9
        h_1(:, :, direction) = - W(direction)*rho_l*h_1_leaf;
    end
% 
    % termo de segunda ordem
    % parte 1
    part_1_second_term = zeros(Nr, Mc); % Nr => y; Mc => x;
    part_1_second_term = diff(rho, 1, 1)*cs2;
    part_1_second_term(Nr, :) = part_1_second_term(Nr - 1, :);
    
    % parte 2
    part_2_second_term = zeros(Nr, Mc); % Nr => y; Mc => x;
    u_r = uy;
    derivada_ux_x = diff(ux, 1, 2);
    derivada_ux_x(:, Mc) = derivada_ux_x(:, Mc - 1);
    derivada_ur_x = diff(u_r, 1, 2);
    derivada_ur_x(:, Mc) = derivada_ur_x(:, Mc - 1);
    derivada_ux_ur_x = derivada_ux_x.*u_r + derivada_ur_x.*ux;
    part_2_second_term = rho_l*derivada_ux_ur_x;

    % parte 3
    part_3_second_term = zeros(Nr, Mc); % Nr => y; Mc => x;
    derivada_ur_r = diff(u_r, 1, 1);
    derivada_ur_r(Nr, :) = derivada_ur_r(Nr - 1, :);
    derivada_ur_ur_r = 2*u_r.*derivada_ur_r;
    part_3_second_term = rho_l*derivada_ur_ur_r;

    % termo parcial para as 9 direcoes
    partial_part_second_term = zeros(Nr, Mc, 9);
    for direction = 1:9
        partial_part_second_term(:,:,direction) = part_1_second_term + ...
        part_2_second_term + part_3_second_term;
    end

    % parte 4
    part_4_second_term = zeros(Nr, Mc, 9); % Nr => y; Mc => x;
    derivada_ux_r = diff(ux, 1, 1);
    derivada_ux_r(Nr, :) = derivada_ux_r(Nr - 1, :);
    for direction = 1:9
        part_4_second_term(:, :, direction) =  ... 
        rho_l*(derivada_ux_r - derivada_ur_x)*C_x(direction);
    end
    
    % parte total
    total_part_second_term = part_4_second_term + partial_part_second_term;

    % calculando matriz de viscosidade
    omega_matrix_second_term = zeros(Nr, Mc, 9);
    viscosity_matrix(1:Nr, 1:Mc) = 3*visc;
    radius = 1:Nr; 
    radius = radius';
    for point_x = 1 : Mc
        viscosity_matrix(:, point_x) = viscosity_matrix(:, point_x)./radius;
    end
    for direction = 1 : 9
        omega_matrix_second_term(:, :, direction) = ... 
        W(direction)*viscosity_matrix;
    end

    % finalmente o termo de segunda ordem
    h_2 = omega_matrix_second_term.*total_part_second_term;

    % colidindo tudo
    f = (1-omega)*f + omega*feq ...
    - sigma_mat9_cima.*(feq- Ft_cima) ...
    - sigma_mat9_esquerda.*(feq- Ft_esquerda) ...
    - sigma_mat9_direito.*(feq- Ft_direito) ...
    - sigma_source.*(feq- Ft_source) ...
    + h_1+ h_2;  

    
    % Ploting the results in real time   
    
    if mod(ta, 100) == 0
%         vorticidade = curl(ux, uy);
%         velocity = sqrt(ux.^2 + uy.^2);
%         imagesc(flip(vorticidade))
        imagesc(flip(rho-1), [-.000001 .000001]);
        axis equal
        pause(.00001);
        disp('Progresso: ');
        disp((ta/total_time*100));
    end
    
    %% Acquisition Data
    
     % Right points
    
    % Plano 5a
    d_plano1 = 565;
    x_5a_r = Mc/2+d_plano1;
    physical_position_4 =( x_5a_r - length_anechoic_condition)*Dx;
    point_distance = radius_pipe/8;    
    for point=1:8
        pressure_r_1(point,ta) = (rho(2+point_distance*(point-1), x_5a_r)-1)*cs2;
        particle_velocity_r_1(point,ta) = ux(2+point_distance*(point-1), x_5a_r);
    end
      
    % Plano 5a + 0.04 m
    d_plano2=565+150;
    x_5a_r2=Mc/2+d_plano2;
    physical_position_5 = (x_5a_r2 - length_anechoic_condition)*Dx;
    for point=1:8
        pressure_r_2(point,ta) = (rho(2+point_distance*(point-1), x_5a_r2)-1)*cs2;
        particle_velocity_r_2(point,ta) = ux(2+point_distance*(point-1), x_5a_r2);
    end
    % Plano 5a + 0.11 m
    d_plano3=565+150+85;
    x_5a_r3=Mc/2+d_plano3;
    physical_position_6 = (x_5a_r3 - length_anechoic_condition)*Dx;
    for point=1:8
        pressure_r_3(point,ta) = (rho(2+point_distance*(point-1), x_5a_r3)-1)*cs2;
        particle_velocity_r_3(point,ta) = ux(2+point_distance*(point-1), x_5a_r3);
    end
    
    % Left points
    
     % Plano 5a
     d_plano4 = -d_plano1;
    x_5a_l1 = Mc/2 + d_plano4;
    physical_position_3 = (x_5a_l1 - length_anechoic_condition)*Dx;
    for point=1:8
        pressure_r_4(point,ta) = (rho(2+point_distance*(point-1), x_5a_l1)-1)*cs2;
        particle_velocity_r_4(point,ta) = ux(2+point_distance*(point-1), x_5a_l1);
    end
    % Plano 5a - 0.04 m
    d_plano5 = -d_plano2;
    x_5a_l2 = Mc/2+d_plano5;
    physical_position_2 = (x_5a_l2 - length_anechoic_condition)*Dx;
    for point=1:8
        pressure_r_5(point,ta) = (rho(2+point_distance*(point-1), x_5a_l2)-1)*cs2;
        particle_velocity_r_5(point,ta) = ux(2+point_distance*(point-1), x_5a_l2);
    end
    % Plano 5a - 0.11 m
    d_plano6 = - d_plano3;
    x_5a_l3 = Mc/2 + d_plano6;
    physical_position_1 = (x_5a_l3 - length_anechoic_condition)*Dx;
    for point=1:8
        pressure_r_6(point,ta) = (rho(2+point_distance*(point-1), x_5a_l3)-1)*cs2;
        particle_velocity_r_6(point,ta) = ux(2+point_distance*(point-1), x_5a_l3);
    end
    
    
    
    
end %  End main time Evolution Loop
%%


% Save pressure in physical units 

 qsi=c_p/cs;
 pressure_1= pressure_r_1*(qsi^2)*1.2;
 pressure_2= pressure_r_2*(qsi^2)*1.2;
 pressure_3= pressure_r_3*(qsi^2)*1.2;
 pressure_4= pressure_r_4*(qsi^2)*1.2;
 pressure_5= pressure_r_5*(qsi^2)*1.2;
 pressure_6= pressure_r_6*(qsi^2)*1.2;
 
 % Saving velocity in physical units 
 particle_velocity_1 = particle_velocity_r_1*qsi;
 particle_velocity_2 = particle_velocity_r_2*qsi;
 particle_velocity_3 = particle_velocity_r_3*qsi;
 particle_velocity_4 = particle_velocity_r_4*qsi;
 particle_velocity_5 = particle_velocity_r_5*qsi;
 particle_velocity_6 = particle_velocity_r_6*qsi;
 
 % organize the velocity and pressure vectors
 
 pressure(1:8,:)=pressure_6;
 pressure(9:16,:)=pressure_5;
 pressure(17:24,:)=pressure_4;
 pressure(25:32,:)=pressure_1;
 pressure(33:40,:)=pressure_2;
 pressure(41:48,:)=pressure_3;
 
 particle_velocity(1:8,:)=particle_velocity_6;
 particle_velocity(9:16,:)=particle_velocity_5;
 particle_velocity(17:24,:)=particle_velocity_4;
 particle_velocity(25:32,:)=particle_velocity_1;
 particle_velocity(33:40,:)=particle_velocity_2;
 particle_velocity(41:48,:)=particle_velocity_3;
 

 % Saving coodirnates
 z = [physical_position_1 physical_position_2 ...
 physical_position_3 physical_position_4 ...
 physical_position_5 physical_position_6];
 r = (((1:8)-1)*point_distance+2)*Dx;
 
 plane_1(1:8,2) = z(1);
 plane_1(1:8,1) = r;
 plane_2(1:8,2) = z(2);
 plane_2(1:8,1) = r;
 plane_3(1:8,2) = z(3);
 plane_3(1:8,1) = r;
 plane_4(1:8,2) = z(4);
 plane_4(1:8,1) = r;
 plane_5(1:8,2) = z(5);
 plane_5(1:8,1) = r;
 plane_6(1:8,2) = z(6);
 plane_6(1:8,1) = r;

 coordinates(1:8,:)=plane_1;
 coordinates(9:16,:)=plane_2;
 coordinates(17:24,:)=plane_3;
 coordinates(25:32,:)=plane_4;
 coordinates(33:40,:)=plane_5;
 coordinates(41:48,:)=plane_6;
 
save stephan_case.mat coordinates pressure particle_velocity