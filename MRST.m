clear
clc
%% input lists
q_oil = []       ; 
q_gas = []       ;
opt_OP = []      ;
opt_rate = []    ;
opt_depth = []   ;
valve_depth = [ 100 , 200 , 300 , 400 , 500 , 600 , 700 , 800 , 900 , 1000 ] ;   % unit = ft
qg_inj = [ 80 , 100 , 200 , 300 , 400 , 500 , 550 , 600 , 700 ] * 1000  ;        % unit = SCF/Day

%% MRST for IPR
for P_i = (0:100)  
    
mrstModule add ad-core ad-blackoil ad-props mrst-gui spe10
gravity reset on

%% Model Geometry
% Grid is a 100-by-1-by-20 Cartesian box with equally sized cells of
% dimensions 25-by-25-by-2.5 feet.
cartDims = [ 100, 1, 20 ];
physDims = cartDims .* [ 25, 25, 2.5 ]*ft;
G = computeGeometry(cartGrid(cartDims, physDims));

%% Petrophysical Properties
% Porosity is constant (=0.2) throughout the formation.  The permeability
% distribution is a correlated geostatistically generated field stored in a
% file supplied by the SPE.
rock = getSPE10_model_1_rock();
clf

%% Define Sources and Sinks
% Model is produced from a single producer located at the far end of the
% model (I==100) constrained at a bottom-hole pressure of 95 Psi.  There is
% a single injector at the near end (I==1) providing pressure support.  The
% injector fills the reservoir with 6.97 cubic metres of gas per day.  Both
% wells have an internal diameter of 1 ft.
IJK = gridLogicalIndices(G);
I   = find(IJK{1}(:,1) == 1);
P1   = find(IJK{1}(:,1) == G.cartDims(1)) ; 
clear IJK

%injection well
W = addWell([], G, rock, I, 'Comp_i', [ 0, 1 ], 'Type', 'rate', ...
            'Val', 6.97*meter^3/day, 'Radius', (1/2)*ft, 'Dir', 'z', ...
            'Sign', +1, 'Name', 'I', 'refDepth', 0*ft);  clear I
   
%vertical well
W = addWell(W, G, rock, P1, 'Comp_i', [ 1, 1 ], 'Type', 'bhp', ...
            'Val', P_i*psia, 'Radius', (1/2)*ft, 'Dir', 'z', ...
            'Sign', -1, 'Name', 'P1', 'refDepth', 0*ft);  clear P1

%% Official Benchmark Relative Permeability Data
% Build a reduced ECLIPSE-style input deck that contains just enough
% information to construct relative permeability curves based on the
% official benchmark data.  In particular we use the fact that the relative
% permeability data is formatted in the same way as ECLIPSE's 'SGOF'
% keyword data.
kr_deck = getSPE10_model_1_relperm();
clf
%% Fluid Properties
% The fluids in this simulation model are incompressible and immiscible
% with constant viscosities.  This means we can use MRST's special purpose
% fluid constructor |initSimpleADIFluid| to create the fluid object.  We
% will however need to use sampled relative permeability curves so we do
% not enter any relative permeability data in this call.
fluid = initSimpleADIFluid('mu'    , [  1, 0.01]*centi*poise, ...
                           'rho'   , [700, 1   ]*kilogram/meter^3, ...
                           'cR'    , 6.0e-5/barsa, ...
                           'phases', 'OG');

%%
% Replace the synthetic relative permeability curves created through
% function |initSimpleADIFluid| with the real benchmark values.
fluid_kr = assignSGOF(fluid, kr_deck.PROPS.SGOF, struct('sat', 1, ...
                                               'interp1d', @interpTable));
fluid.krG = fluid_kr.krG{1};
fluid.krO = fluid_kr.krOG{1};              clear fluid_kr

%% Form Reservoir Model
% This is an incompressible, immiscible oil/gas system.
model = GenericBlackOilModel(G, rock, fluid, 'gravity', gravity, 'disgas', false,...
    'vapoil', false, 'water', false, 'oil', true, 'gas', true);

%% Initialise Formation
% Formation is initially filled with oil and the initial pressure at the
% top of the model is 100 Psi.
region = getInitializationRegionsBlackOil(model, 0, 'datum_pressure', 100*psia);
state0 = initStateBlackOilAD(model, region);
clf
%%
timesteps = [ 0.1, 0.2 ]*day;

% Set up the schedule containing both the wells and the timesteps
schedule = simpleSchedule(timesteps, 'W', W);

[wellSols, states] = ...
   simulateScheduleAD(state0, model, schedule);

% unit = m^3/ day
q_oil = [ q_oil , -wellSols{1}(2).qOs] ;  
q_gas = [ q_gas , -wellSols{1}(2).qGs] ;

end

P_avg =  mean(state0.pressure) / 6895 ;                                    % average pressure  unit = psi
q_oil = q_oil * 3600 *24 / ( 0.3048^3 * 5.615)  ;                          % m^3/d to STB/d
q_gas = q_gas * 3600 *24 / ( 0.3048^3 * 5.615) * 5.615 ;                   % m^3/d to SCF/d
P_wf = (0:100)  ;

%% Gas Lift :)
dp_min_valve = 70   ;            % minimum differential pressure for valve   unit = psi
dp_max_valve = 3000 ;            % maximum differential pressure for valve   unit = psi
HHP = 40 ;                       % Hydraulic horsepower                      unit = hp
tubing_dia = 10   ;              % Tubing diameter                           unit = inch
casing_dia = 30   ;              % Casing diameter                           unit = inch
gas_vis = 0.02    ;              % gas viscosity                             unit = cp
rho_g = 0.042     ;              % gas density                               unit = lb / ft^3
Rough = 0.0006    ;              % Pipe Roughness                            unit = in
eps = Rough / casing_dia   ;     % epsilon / d
P_inlet_comp = 60   ;            % Compressor inlet pressure                 unit = psi
P_inj = P_inlet_comp * ( 10^4 * HHP ./ ( 2.23 * qg_inj ) + 1 ).^5 ;          % injection pressure   unit = psi
for i = ( 1 : length(valve_depth) )
    for j = ( 1 : length(qg_inj) )     
        K_valve = 8.5   ;
        valve_dia = 1.5   ;                                                % Valve diameter       unit = in
        delta_p_energy = (K_valve / (2 * valve_dia^4)) * rho_g * (qg_inj(j))^2 * 6 * 10^(-10) ; % unit = psi 
        v_inj = ( 144 * 4 * qg_inj(j) ) / ( pi * 86400 * casing_dia^2 )   ;                     % unit = ft / s
        Re = 1488 * rho_g * v_inj * casing_dia / ( 12 * gas_vis )     ;    % Reynolds number
        if Re < 2000
            f = 64 / Re   ;
        else
            f = 1 / ( -2 * log10 ( eps / 3.7 - ( 5.02 / Re ) * log10( eps / 3.7 + 13 / Re )))^2 ; % Zigrang , Sylvester
        end
        
        dp = (- rho_g * (v_inj)^2 * f / ( 2 * casing_dia / 12 ) + rho_g ) / 144 ; %Pressure gradient unit=psi/ft
        P_valve_u = dp * valve_depth(i) + P_inj(j)   ;                            
        
        % CPR
        P_wh = [] ;           % Well head Pressure list                    unit = psi
        S = 32   ;            % choke size                                 unit = 1/64 in      
        for k = ( 1 : length(q_oil) ) 
            q_gas(k) =  2 * qg_inj(j)   ;                                  % total gas  unit = SCF / Day
            R = q_gas(k) / q_oil(k)               ;                        % GOR        unit = SCF / STB
            P_1 = 10 * ( R ^ 0.546 ) * q_oil(k) / S ^ 1.89   ;
            P_wh = [ P_wh , P_1 ]  ;
        end
        
        % TPR  by using Beggs and brill

        % TPR #1 : Wellhead -----> Valve
        
        h = valve_depth(i) ;              % Tubing shoe depth              unit = ft
        WC = 0   ;                        % Water cut
        t = 1  ;
        P_wf_beggs = []  ;                % P_wf by beggs and brill        unit = Psi
        while t <= length(P_wh) 

        q_L = q_oil(t)   ;                % Liquid production rate         unit = stb/day
        GLR = q_gas(t) / q_oil(t)   ;     % Gas liquid ratio               unit = scf / stb
        P_head = P_wh(t)      ;           % Well head pressure             unit = psi
        P_valve_d = real (poettmann_carpenter (P_head,tubing_dia,h,q_L,WC,GLR)) ;
        delta_p_valve = P_valve_u - P_valve_d   ;
        
        if (dp_min_valve < delta_p_valve) && (delta_p_valve < dp_max_valve)
            % TPR #2 : Valve -----> bottom hole
            h = 1050 - valve_depth(i) ;          % Tubing shoe depth       unit = ft
            P_head = P_valve_d      ;            % Well head pressure      unit = psi
            GLR = 0   ;                          % Gas liquid ratio        unit = scf / stb
            P_wf_b = real(poettmann_carpenter (P_head,tubing_dia,h,q_L,WC,GLR)) ;
            P_wf_beggs = [ P_wf_beggs , P_wf_b ]   ;
            t = t + 1   ;
            
        else
            P_wf_beggs = [ -999 ]   ;  
            t = length(P_wh) + 2   ;
        end

        end
        
        if length(P_wf_beggs) == 1
            P_wf_beggs = []   ;
        else
            Dif = P_wf - P_wf_beggs    ;
            for u = (1:length(Dif)-1)
                if Dif(u) == 0
                    OP1 = q1(u)  ;
                elseif Dif(u) * Dif(u+1) < 0
                    OP1 = ((abs(Dif(u+1)))/(abs(Dif(u)) + abs(Dif(u+1)))) * q_oil(u) + ((abs(Dif(u)))/(abs(Dif(u)) + abs(Dif(u+1)))) * q_oil(u+1) ;
                end
            end

        end
        
        
         
        if  exist('OP1','var') == 1
            opt_OP = [ opt_OP , OP1 ]   ;
            opt_rate = [ opt_rate , qg_inj(j) ]   ;
            opt_depth = [ opt_depth , valve_depth(i) ]   ; 
            mm = i + ( j - 1 ) * length(valve_depth)   ;
            figure(mm+1)
            plot(q_oil,P_wf_beggs,'b') ;
            hold on
            plot( q_oil , P_wf,'k') ;
            title(['Nodal Analysis || Valve depth = ',num2str(valve_depth(i)),' - Injection rate = ',num2str(2*qg_inj(j))]) ; 
            legend('TPR','IPR')  ;
            ylabel('P_wf (Psi)')  ;
            xlabel('Flow rate (bbl/day)')  ;
            
           
        end
         
        clear OP1
        clear Dif
    end
end

max = -1   ;
ind = -1   ;
for i = (1: length(opt_OP))
    if opt_OP(i) > max
        max = opt_OP(i)   ;
        ind = i  ;
    end 
end

disp(['Operating Point = ',num2str(max),' bbl/day'])
disp(['Optimum injection rate = ',num2str(2 * opt_rate(ind)),' Scf/day'])
disp(['Optimum injection point = ',num2str(opt_depth(ind)),' ft'])
 
%% TPR function

function [ bhp ] = poettmann_carpenter (P_head,D,h,q_L,WC,GLR) 
%inputs:
% P_head = Well head pressure              unit = psi
% D = Tubing diameter                      unit = in
% h = Tubing shoe depth                    unit = ft
% q_L = Liquid production rate             unit = stb/day
% WC = Water cut
% GLR = Gas liquid ratio                   unit = scf/stb

T_head = 100     ;             % T_head = Well head temperature            unit = F
T_bh = 200       ;             % T_bh = Bottom hole temperature            unit = F
API_o = 70.64    ;             % API_o = Oil gravity                       unit = API
gamma_w = 1.05   ;             % gamma_w = Water specific gravity 
gamma_g = 0.62   ;             % gamma_g = Gas specific gravity
rho_air = 0.076  ;             % air density                               unit = lb/ft^3
Bw = 1.2         ;             % Formation volume factor for water         unit = rb / stb

% Calculation
gamma_o = 141.5 / (131.5 + API_o )   ;                                     % Oil specific gravity
D = D / 12   ;                                                             % in to ft
q_o = ( 1 - WC ) * q_L  ;      % Oil production rate                       unit = stb/day
q_w = q_L - q_o   ;            % Water production rate                     unit = stb/day
WOR = q_w / q_o   ;            % Water oil ratio
GOR = GLR / ( q_o / (q_o + q_w ))   ;    % Gas oil ratio                   unit = scf / stb
% Mass associated with 1 stb of oil 
M = 350.17 * ( gamma_o + WOR * gamma_w ) + GOR * gamma_g * rho_air   ;     % unit = lb
D_rho_v = ( 1.4737 * 10^(-5) * M * q_o ) / D   ;      % Inertial force     unit = lb / day-ft
f_2F = 4 * 10^(1.444 - 2.5 * log10(D_rho_v))   ;                           % Two-phase friction factor
K_ave = ( f_2F * (q_o * M )^2 ) / ( 7.4137 * 10 ^10 * D^5 )   ;

% Well head :
% Solution gas oil ratio at well head             unit = scf/stb
R_s = 0   ;
% Oil formation volume factor at wellhead         unit = rb/stb
B_o = 0.9759 + 0.000012 * ( R_s * sqrt ( gamma_g / gamma_o ) + 1.25 * T_head )^1.2           ; 

%Z factor Calculation at wellhead
Ppc = 756.8 - 131.07  *gamma_g - 3.6 * gamma_g^2     ;                     % Pseudocritical Pressure
Tpc = 169.2 + 349.5* gamma_g - 74 * gamma_g^2        ;                     % Pseudocritical Temperature
Pr = P_head / Ppc                                    ;                     % Pseudoreduced Pressure
Tr = ( T_head + 460 )/ Tpc                           ;                     % Pseudoreduced Temperature
A1 = 1.39 * ( Tr - 0.92 )^0.5 - 0.36 *Tr - 0.101     ; 
B1 = ( 0.62 - 0.23 * Tr ) * Pr + ( 0.066 / ( Tr - 0.86 ) - 0.037 ) * Pr^2 + 0.32 / 10^( 9 * ( Tr - 1 )) * Pr^6   ; 
C1 = 0.132 - 0.32 * log10( Tr )                      ;  
D1 = 10^( 0.3106 - 0.49 * Tr + 0.1824 * ( Tr^2 ) )   ; 
Z_head = A1 + ( 1 - A1 ) / exp( B1 ) + C1 * Pr^D1    ;                     % Z Factor

% Volume associated with 1 stb of oil at wellhead 
V_m = 5.615 * (B_o + WOR * Bw ) + ( GOR - R_s ) * ( 14.7 * Z_head * ( T_head + 460 )) / ( P_head * 520 )   ;
rho1 = M / V_m   ;                  % Fluid density at wellhead            unit = lb / ft^3

% Bottom hole :
PD = P_head + ( rho1 * h * 144 )    ;      %Estimated Pressure Drop         unit = psi
P_w = 0   ;
while abs ( PD - P_w ) > 0.01  
    P_w = PD   ;
    % Solution gas-oil ratio at bottom hole          unit = scf/stb
    R_s = 0   ;
    % Oil formation volume factor at bottom hole     unit = rb / stb
    B_o = 0.9759 + 0.000012 * ( R_s * sqrt ( gamma_g / gamma_o ) + 1.25 * T_bh )^1.2       ;  
    %Z factor Calculation at bottom hole
    Ppc = 756.8 - 131.07  *gamma_g - 3.6 * gamma_g^2     ;                     % Pseudocritical Pressure
    Tpc = 169.2 + 349.5* gamma_g - 74 * gamma_g^2        ;                     % Pseudocritical Temperature
    Pr = PD / Ppc                                    ;                         % Pseudoreduced Pressure
    Tr = ( T_bh + 460 )/ Tpc                           ;                       % Pseudoreduced Temperature
    A1 = 1.39 * ( Tr - 0.92 )^0.5 - 0.36 *Tr - 0.101     ; 
    B1 = ( 0.62 - 0.23 * Tr ) * Pr + ( 0.066 / ( Tr - 0.86 ) - 0.037 ) * Pr^2 + 0.32 / 10^( 9 * ( Tr - 1 )) * Pr^6   ; 
    C1 = 0.132 - 0.32 * log10( Tr )                      ;  
    D1 = 10^( 0.3106 - 0.49 * Tr + 0.1824 * ( Tr^2 ) )   ; 
    Z_bh = A1 + ( 1 - A1 ) / exp( B1 ) + C1 * Pr^D1    ;                       % Z Factor
    % Volume associated with 1 stb of oil at bottom hole
    V_m = 5.615 * (B_o + WOR * Bw ) + ( GOR - R_s ) * ( 14.7 * Z_bh * ( T_bh + 460 )) / ( PD * 520 )   ;
    rho2 = M / V_m   ;                   % Fluid density at bottom hole    unit = lb / ft^3
    rho_ave = ( rho1 + rho2 ) / 2    ;   % The average fluid density       unit = lb / ft^3
    delta_P = ( rho_ave + K_ave / rho_ave ) * h / 144   ;                % unit = psi
    PD = delta_P + P_head   ;                                          
end

bhp = PD   ;

end