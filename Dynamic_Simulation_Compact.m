%function [score] = Dynamic_Simulation_Compact(p,e,k,k_reduced,reduced_index,feed_start,def_C_fit)
  % ,V_All, p, Conc_met, R, C_fit, norm_C_fit, linear_C_fit, time_vector
%INPUTS
%p: process variables 
%e: genetic variables
%c: cell line where c = 1 for cell line A; c = 2 for cell line B
% OUTPUTS
% A: matrix of simulated time, fluxes, and metabolite concentrations

%% FOR PERFORMANCE PROFILING - Kevin
% path_from_base='../';
% default_values_path='mat_files/nolan_defaults.mat'; 
% c_fit_path='mat_files/def_c_fit.mat'; 
% loader=load([path_from_base,c_fit_path]); 
% def_C_fit=loader.def_C_fit;
%k(reduced_index)=k_reduced; 
load('nolan_defaults.mat')
p=default_nolan.p;
e=default_nolan.e;
k=default_nolan.k;
%% Initialization

debugger_original_container=struct;
%     if c == 1
%         %k = [3.8;2;1;3;1;5;2.2;6;3.5;1.5;4;2.5;5;2;2.2;2;2.5;1;1;5;1;1;7;4.75;2;3;1;2;1.5;5.5;0.9;2.5;1;3;3;1;2.5;4;3;2;3;5.25;3;5;1;2;2;];
%         k;
% % A parameters
%     end
%     if c == 2
%         k = [3.8;3;1;3;2;5;2.3;4;6.25;2;5;1.8;2.5;2;1.5;3;4;3;1;5;1;1;7;8.8;2;5;2;2;1.5;6.5;1;2.5;1;3;2.7;1;1.8;5;3;2;2.2;4.5;4;9;0.4;1;2;];
% B parameters
%     end

    c=1;
    R = 6.2; %Redox variable to accomodate the effect of NADH on reaction rate kinetics

    i = 1;
    t_interval = 4;
    t = (i-1)/t_interval; %time interval is 1/4th of a day
    V = [];
    initial_temp = 37;
%------------------------------------------------------------------------
% process variables
% A default: [31,3,1.80,10,0.03];
% B default: [31,3,1.40,10,0.018];
    shift_temp = p(1); % the T the culture is reduced to or shifted to at a specified time following inoculation; default 31
    shift_day = p(2); % the day the temperature is reduced ;default 3
    seed_density = p(3); % cell density at time of inoculation; default 1.80
    harvest_day = p(4); % the day the culture is terminated ;default 10
    feed_daily = p(5); % default 0.03 (i.e. 3% v/v)
    feed_start = 3; %feed_daily starts on day 3
% genetic variables
% -1 = down-regualtion, reset to 1/scale factor
% 0 = no genetic modulation, reset to 1
% 1 = up-regulation, reset to scale factor
%k is matrix of the 47 kinetic parameters explained from line ~129-179.
%Based on genetic regulation of each reaction, parameters are adjusted
    k(1) = k(1)*e(1);
    k(7) = k(7)*e(2);
    k(9) = k(9)*e(3);
    k(10) = k(10)*e(3);
    k(15) = k(15)*e(4);
    k(16) = k(16)*e(5);
    k(21) = k(21)*e(6);
    k(22) = k(22)*e(6);
    k(24) = k(24)*e(7);
    k(25) = k(25)*e(7);
    k(30) = k(30)*e(8);
    k(31) = k(31)*e(9); % vmax12f
    k(32) = k(32)*e(9); % vmax12r
    k(35) = k(35)*e(10); % vmax13
    if c == 1
        initial_BIOM = seed_density*2.31;
        initial_ANTI = 0.14;
        initial_GLC = 72.61;
        initial_LAC = 5.55;
        initial_ALA = 0.67;
        initial_ASN = 20.27;
        initial_ASP = 2.05;
        initial_CYS = 1.50;
        initial_GLN = 3.90;
        initial_GLY = 3.30;
        initial_SER = 11.06;
        initial_NH3 = 0.68;
        initial_GLU = 0.62;
        initial_VCD = seed_density;
        initial_TCD = seed_density;
        initial_DCD = 0.00;
    end
    if c == 2
        initial_BIOM = seed_density*2.31;
        initial_ANTI = 0.06;
        initial_GLC = 70.00;
        initial_LAC = 5.11;
        initial_ALA = 0.62;
        initial_ASN = 18.87;
        initial_ASP = 1.82;
        initial_CYS = 1.50;
        initial_GLN = 4.00;
        initial_GLY = 3.15;
        initial_SER = 10.32;
        initial_NH3 = 0.52;
        initial_GLU = 0.47;
        initial_VCD = seed_density;
        initial_TCD = seed_density;
        initial_DCD = 0.00;
    end
    initial_concen = [initial_BIOM; initial_ANTI; initial_GLC; initial_LAC;initial_ALA; initial_ASN; initial_ASP; initial_CYS; initial_GLN;initial_GLY; initial_SER; initial_NH3; initial_GLU; initial_VCD;initial_TCD; initial_DCD];
    C = initial_concen;
    feed_BIOM = 0.00;
    feed_ANTI = 0.00;
    feed_GLC = 444.4;
    feed_LAC = 0.00;
    feed_ALA = 6.00;
    feed_ASN = 54.01;
    feed_ASP = 15.04;
    feed_CYS = 4.68;
    feed_GLN = 0.00;
    feed_GLY = 6.00;
    feed_SER = 45.16;
    feed_NH3 = 0.00;
    feed_GLU = 6.00;
    feed_VCD = 0.00;
    feed_TCD = 0.00;

    feed_concen = [feed_BIOM; feed_ANTI; feed_GLC; feed_LAC; feed_ALA;feed_ASN; feed_ASP; feed_CYS; feed_GLN; feed_GLY; feed_SER; feed_NH3;feed_GLU; feed_VCD; feed_TCD];
    feed_day = feed_start:harvest_day;
    feed_percent(1:length(feed_day)) = feed_daily;
%------------------------------------------------------------------------
    %scaling back parameter values to their original OOM (OOM of parameters
    %was changed for initial optimization of parameters, not in this code)
    vmax1 = 1000*k(1);
    Ki1 = 10*k(2);
    Km1 = 10*k(3);
    exp1a = k(4);
    exp1b = k(5);
    TC1b = k(6);
    vmax2 = 1000*k(7);
    Km2 = k(8);
    vmax3f = 100*k(9);
    vmax3r = 100*k(10);
    Km3a = k(11);
    Km3b = k(12);
    Km3c = k(13);
    TC3 = k(14);
    vmax8f = 1000*k(15);
    vmax8r = 100*k(16);
    Km8a = k(17);
    Km8b = k(18);
    Km8c = k(19);
    TC8b = k(20);
    vmax9f = k(21);
    vmax9r = k(22);
    Km9 = 0.1*k(23);
    vmax10f = 100*k(24);
    vmax10r = 10*k(25);
    Km10a = 0.1*k(26);
    Km10b = k(27);
    Km10c = k(28);
    TC10b = k(29);
    vmax11 = 0.1*k(30);
    vmax12f = 0.1*k(31);
    vmax12r = 10*k(32);
    Km12a = k(33);
    Km12b = k(34);
    vmax13 = 10*k(35);
    Km13 = k(36);
    vmax16 = 1000*k(37);
    Km16a = 0.1*k(38);
    Km16b = 0.1*k(39);
    Km16c = 10*k(40);
    TC16b = k(41);
    vmax17 = 100*k(42);
    Ki17 = 10*k(43);
    exp17a = 0.1*k(44);
    exp17b = k(45);
    vmax33a = 0.10*k(46);
    vmax33b = 0.10*k(47);
    VCD_peak = 1;
    LAC_peak = 1;
    R_peak = 0.5;
    q_stop = 0;
%     TC2 = 31;

    %% Main loop
while i <= (t_interval*harvest_day+1) %i goes from 1 to 41
    BIOM = C(1,:);
    GLC = C(3,:);
    LAC = C(4,:);
    ALA = C(5,:);
    ASN = C(6,:);
    ASP = C(7,:);
    CYS = C(8,:);
    GLN = C(9,:);
    GLY = C(10,:);
    SER = C(11,:);
    NH3 = C(12,:);
    GLU = C(13,:);
    VCD = C(14,:);

%% Temp Calculations
%     disp(['iteration: ',num2str(i),'. ',datestr(clock,21)]);
    temp = initial_temp - ((initial_temp - shift_temp)/(1 + 1*exp(-10*(t(i) - (shift_day + 0.25)))));
    leak = 1 + (1.5/(1 + 1*exp(-10*(t(i) - (shift_day + 0.25)))));
    Ta = (0 - 1)*(initial_temp - temp)/(initial_temp - shift_temp) + 1;
    Tb = 1 - Ta;
%--------------------------------------------------------------------
    TC1bx = Arrhenius(TC1b,shift_temp,1);
    TC1 = (TC1bx - 1)*(initial_temp - temp)/(initial_temp - shift_temp) + 1;
%     TC2x = Arrhenius(TC2,shift_temp,1);
    TC3x = Arrhenius(TC3,shift_temp,1);
    TC8bx = Arrhenius(TC8b,shift_temp,1);
    TC8 = (TC8bx - 1)*(initial_temp - temp)/(initial_temp - shift_temp) + 1;
    TC10bx = Arrhenius(TC10b,shift_temp,1);
    TC10 = (TC10bx - 1)*(initial_temp - temp)/(initial_temp -shift_temp) + 1;
    TC16bx = Arrhenius(TC16b,shift_temp,1);
    TC16 = (TC16bx - 1)*(initial_temp - temp)/(initial_temp -shift_temp) + 1;



    exp1bx = ArrheniusExp(exp1b,shift_temp,exp1a);
    exp1 = (exp1bx - exp1a)*(initial_temp - temp)/(initial_temp -shift_temp) + exp1a;
    exp17bx = ArrheniusExp(exp17b,shift_temp,exp17a);
    exp17 = (exp17bx - exp17a)*(initial_temp - temp)/(initial_temp -shift_temp) + exp17a;

    capture_variables={'temp','leak','TC1bx','TC1','TC3x','TC8bx','TC8',...
        'TC10bx','TC10','TC16bx','TC16','exp1bx','exp1','exp17bx','exp17'};
    for j=1:length(capture_variables)
        debugger_original_container.(capture_variables{j})(i)=eval(capture_variables{j});
    end

    %% Kinetic Flux Calculations for intracellular reactions
    V(1,i) = ((vmax1/TC1)*(1/(1+(LAC(i)/Ki1)^exp1))*(GLC(i)/Km1))/(1 +(GLC(i)/Km1));
%--------------------------------------------------------------------
    if R(i) >= 1
        V(2,i) = (R(i) - 1)*vmax2;
%         {i, 'eqn_1'}
%         V(2,i);
    else
%         {i,'eqn_2'}
        V(2,i) = ((R(i) - 1)*vmax2*(LAC(i)/Km2))/(1 + (LAC(i)/Km2));
    end
%--------------------------------------------------------------------

    if t(i) <= shift_day
        V(3,i) = (vmax3f*(GLC(i)/Km3c) - vmax3r*(ALA(i)/Km3b))/(1 +(GLC(i)/Km3c) + (ALA(i)/Km3b));
    else
    V(3,i) = (vmax3f*(LAC(i)/Km3a) - vmax3r*(ALA(i)/Km3b))/(1 +(LAC(i)/Km3a) + (ALA(i)/Km3b));
    end
    if (t(i) > 1) && (V(3,i-1)<0)
        V(3,i) = -(vmax3r/TC3x)*(ALA(i)/Km3b)/(1 + (ALA(i)/Km3b));
    end
%--------------------------------------------------------------------

    V(8,i) = ( (vmax8f/TC8)*(GLN(i)/Km8a) - vmax8r*(GLU(i)/Km8b) *(NH3(i)/Km8c) )/...
        ( 1 + (GLN(i)/Km8a) + (GLU(i)/Km8b) +(NH3(i)/Km8c) + (GLU(i)/Km8b)*(NH3(i)/Km8c) );
%--------------------------------------------------------------------

    V(10,i) = ((vmax10f/TC10)*(ASN(i)/Km10a) - (vmax10r/TC10) *(ASP(i)/Km10b)*(NH3(i)/Km10c))/ ...
        ...
        (1 + (ASN(i)/Km10a) +(ASP(i)/Km10b)+ (NH3(i)/Km10c) + (ASP(i)/Km10b)*(NH3(i)/Km10c));
%--------------------------------------------------------------------

    V(11,i) = vmax11*V(10,i);
%--------------------------------------------------------------------

    V(9,i) = (vmax9f*V(3,i)*(NH3(i)/Km9) - vmax9r*V(11,i))/(1+(NH3(i)/Km9));
%--------------------------------------------------------------------

    V(13,i) = Ta*((R(i)*vmax13*(CYS(i)/Km13))/(1 + (CYS(i)/Km13))) +Tb*((vmax13*(CYS(i)/Km13))/(1 + (CYS(i)/Km13)));
%--------------------------------------------------------------------

    V(16,i) = Ta*(((vmax16/TC16)*(GLN(i)/Km16a)*(ASN(i)/Km16b))/(1 +(GLN(i)/Km16a) + (ASN(i)/Km16b) + (GLN(i)/Km16a)*(ASN(i)/Km16b))) + Tb*(((vmax16/TC16)*(GLC(i)/Km16c)*(ASN(i)/Km16b))/(1 + (GLC(i)/Km16c) + (ASN(i)/Km16b) + (GLC(i)/Km16c)*(ASN(i)/Km16b)));
%--------------------------------------------------------------------

    V(12,i) = (vmax12f*V(16,i)*(SER(i)/Km12a) - vmax12r*((GLY(i)/Km12b)^2))/(1 + (SER(i)/Km12a) + (GLY(i)/Km12b) + (GLY(i)/Km12b)^2);
%--------------------------------------------------------------------
    V(17,i) = vmax17/(1 + (LAC(i)/Ki17)^exp17);
%--------------------------------------------------------------------
    V(33,i) = vmax33a*V(8,i) + vmax33b*V(11,i);
%--------------------------------------------------------------------
    V(18,i) = V(16,i);
    V(19,i) = V(17,i);
    V(22,i) = V(3,i) - (0.0837953950666345*V(16,i) +0.061377245508982*V(17,i));
%--------------------------------------------------------------------

%% Stoichiometric flux FBA + minimization
    vmeas_index = [1,2,3,8,9,10,11,12,13,16,17,18,19,33]';
    vmeas = V(vmeas_index,i);
    debugger_original_container.vmeas(:,i)=vmeas;
    err = 0.10;
    [v,debugger_original_container.ub(:,i),debugger_original_container.lb(:,i)] = FBA(vmeas_index,vmeas,t(i),shift_day,leak,err);
    q = 1; %q is the iteration counter for FBA code

    %% Not used in default
    if q_stop < 1
        q_stop;
        while ((v(1) > 10000) || (v(2) < -10000) || (v(8) < -5000) ||(v(17) < 150)) && (q < 31)
%             {'Conditional',i,q_stop,q}
            v = FBA(vmeas_index,vmeas,t(i),shift_day,leak,err);
            q = q + 1;
            if q >= 10
                err = 1; %initial "error" was 0.1 but if code doesn't solve in 10 iterations, change to 1 to enhance range of lb and ub (see FBA script)
            end
            if q > 30
                q_stop = 1; %if code doesn't solve in 30 ietrations, stop
            end
        end
    else
%         {'Normal',i}
        v = FBA(vmeas_index,vmeas,t(i),shift_day,leak,err);
    end


    V_All(:,i) = v;
    %v(2)
%     V_FBA(:,i) = v(vmeas_index);
%     V_Ext(:,i) = v([18:29,33]);

    %% Time step
    i = i+1;
    t(i) = (i-1)/t_interval;
%     t(i) = TIME(i);

    %% Update redox / concentrations

    R(i) = (2*v(1) + 0.6391*v(16))/v(34); % Redox

    C(1,i) = C(1,i-1)+v(18)*(t(i)-t(i-1))*C(14,i-1)/1000; %BIOM, +rxn18

    C(2,i) = C(2,i-1)+v(19)*(t(i)-t(i-1))*C(14,i-1)/1000; %ANTI, +rxn19

    C(3,i) = C(3,i-1)-v(20)*(t(i)-t(i-1))*C(14,i-1)/1000; %GLC, -rxn20

    C(4,i) = C(4,i-1)+v(21)*(t(i)-t(i-1))*C(14,i-1)/1000; %LAC, +rxn21

    C(5,i) = C(5,i-1) + v(22)*(t(i)-t(i-1))*C(14,i-1)/1000; %ALA, +rxn22

    C(6,i) = C(6,i-1) - v(23)*(t(i)-t(i-1))*C(14,i-1)/1000; %ASN -rxn23

    C(7,i) = C(7,i-1) + v(24)*(t(i)-t(i-1))*C(14,i-1)/1000; %ASP, +rxn24

    C(8,i) = C(8,i-1) - v(25)*(t(i)-t(i-1))*C(14,i-1)/1000; %CYS, -rxn25

    C(9,i) = C(9,i-1) - v(26)*(t(i)-t(i-1))*C(14,i-1)/1000; %GLN, -rxn26

    C(10,i) = C(10,i-1) + v(27)*(t(i)-t(i-1))*C(14,i-1)/1000; %GLY, +rxn27

    C(11,i) = C(11,i-1) - v(28)*(t(i)-t(i-1))*C(14,i-1)/1000; %SER, -rxn28

    C(12,i) = C(12,i-1) + v(29)*(t(i)-t(i-1))*C(14,i-1)/1000; %NH3, +rxn29

    C(13,i) = C(13,i-1) + (v(25) + v(33))*(t(i)-t(i-1))*C(14,i-1)/1000; %GLU, +rxn25 +rxn33

    C(15,i) = C(1,i)/2.31; %overall cell density


    %% Viability and VCD calculations
    if c == 1
        via = 1.07 - (0.80/((1 + exp(-0.35*(t(i) - 2.50)))^(1/0.25)));% A
    end
    if c == 2
        via = 1.10 - (0.90/((1 + exp(-0.38*(t(i) - 2.50)))^(1/0.25)));% B
    end
    debugger_original_container.viability(i-1)=via;
    if via > 1
        via = 1;
    end
%     t(i)


    C(14,i) = C(15,i)*via; %VCD = viability * overall density
    C(16,i) = C(15,i) - C(14,i); 

% Feeding, avoiding negative concentrations, setting variable peaks
    if C(4,i) > LAC_peak
        LAC_peak = C(4,i);
    end
    if C(14,i) > VCD_peak
        VCD_peak = C(14,i);
    end
    if (R(i) > R_peak) && (t(i) > 5)
        R_peak = R(i);
    end
    if isempty(find(C(:,i)<0, 1)) < 1 % adjust if concentration is <=0
%         {i,'All',find(C(:,i)<0)}

        C(C(:,i)<0,i) = 0.01;
    end
    if C(3,i) < 1 % add extra glucose if concentration is too low
%         {i,'glc'}
        C(3,i) = 5;
    end
    if (C(5,i) < 0.25) && (t(i) > 1) % add alanine if concen. is too low
%         {i,'ala'}
        C(5,i) = 1;
    end
    if C(6,i) < 0.25 % add extra asparagine if concentration is too low
%         {i,'asn'}
%         C(6,i)
        C(6,i) = 1;
    end
    if C(7,i) < 0.25 % add extra aspartate if concentration is too low
%         {i,'asp'}
        C(7,i) = 1;
    end
    if C(8,i) < 0.1 % add extra cystine if concentration is too low
%         {i,'cys'}
%         {C(8,i),i}
        C(8,i) = 0.2;
    end
    if C(10,i) < 0.25 % add extra glycine if concentration is too low
%         {i,'gly'}
        C(10,i) = 1;
    end
    if C(11,i) < 0.25 % add extra serine if concentration is too low
%         {i,'ser'}
        C(11,i) = 1;
    end
    if isempty(intersect(t(i),feed_day)) < 1 % 0 if feed day
%         {i,t(i)}
        percent_index = find(feed_day==t(i));
        C(1:15,i) = C(1:15,i)*(1 - feed_percent(percent_index)) + feed_concen*feed_percent(percent_index);

    end
end

%% Output
A = [t(1:t_interval*harvest_day+1); V_All; C(:,1:t_interval*harvest_day+1)];

Conc_met = C(:,1:t_interval*harvest_day+1);

C_fit = [Conc_met(1:14,:);R(1:length(Conc_met(1,:)))];

ANTI_final = C(2,t_interval*harvest_day);

LAC_final = C(4,t_interval*harvest_day);

vLAC_final = sign(V(2,t_interval*harvest_day));

LAC_final = LAC_final*vLAC_final;

NH3_final = C(12,t_interval*harvest_day);

R_final = R(t_interval*harvest_day);

IVCD_final = sum(C(14,:))*(1/t_interval);

Qp_final = ANTI_final/IVCD_final;

Z = [ANTI_final, LAC_final, NH3_final, R_final, IVCD_final, Qp_final, LAC_peak, VCD_peak, R_peak, via];

%% Added for fitting code
counter=1;
linear_C_fit=nan(numel(C_fit),1);

% Normalize metabolite concentration by metabolite's max conc
norm_C_fit=C_fit';
C_fit=C_fit';

for j=1:size(norm_C_fit,2)
    norm_C_fit(:,j)=norm_C_fit(:,j)/max(norm_C_fit(:,j));
end

% Linearize vector
% normalize_output=1;
% if normalize_output
%     C_fit_reporter=norm_C_fit;
% else
%     C_fit_reporter=C_fit';
% end
% for j = 1:size(C_fit_reporter,2)
%     for i = 1:size(C_fit_reporter,1)
%         linear_C_fit(counter) = C_fit_reporter(i,j);
%         counter = counter + 1;
%     end
% end


time_vector=[1:t_interval*harvest_day+1]';

%score=sum(sum((abs(C_fit-def_C_fit)./def_C_fit).^2));

%end

%% Helper functions

function TC_new = Arrhenius(TC_old,shift_temp,start)
    P = polyfit([1/310 1/304],[log(start) log(1/TC_old)],1);
    k = exp(P(2))*exp(P(1)./(shift_temp + 273));
    TC_new = 1./k;
end

function TC_new = ArrheniusExp(TC_old,shift_temp,start)
    P = polyfit([1/310 1/304],[log(start) log(TC_old)],1);
    k = exp(P(2))*exp(P(1)./(shift_temp + 273));
    TC_new = k;
end

function [x,ub,lb] = FBA(~,vmeas,~,~,leak,err)
    warning off;
    vmeas(2);
    vmeas_index = [1,2,3,8,9,10,11,12,13,16,17,18,19,33]';

%     rand('state',sum(100*clock));
rng('default');
    A = [];
    b = [];
    % Stoichiometric matrix setup
%     Aeq = [0,0,1,1,-1,0,0,0,-1,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,-0.0838,-0.0614,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,-0.041,-0.0344,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,1,-1,0,0,0,0,-0.0804,-0.0389,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0;3,0,0,0,1,0,0,0,0,0,0,0,0,2.5,1.5,-8.6825,-9.2,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,2,1,0,1,0,0,0,0,-1,0,0,0,1.509,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,-0.0261,-0.024,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2,0,0;0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0;0,0,0,0,1,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.50;-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.452,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,-1,0,0,0,0,1,1,0,1,0,0,0,0,0.0082,-0.0479,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,-1,0;0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,-0.0873,-0.0449,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,-0.0807,-0.0719,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0;0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,1,-1,-1,0,0,0,0,0,0,0,0,0.4445,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;2,-1,0,0,0,0,0,0,0,0,0,-1,-1,0,0,0.6391,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1;0,0,0,2,1,1,0,0,-1,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.50;0,0,0,0,0,0,0,1,-1,1,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,-0.50,-0.50,-0.1002,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,-1,0,0;0,0,0,-1,0,1,0,0,0,0,1,0,0,0,0,-0.4270,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;2,-1,-1,-1,0,0,1,0,0,0,0,0,0,0,0,0.2085,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,-0.07130,-0.1000,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0;];
    Aeq = [0,0,1,1,-1,0,0,0,-1,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,-0.0838,-0.0614,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0;...
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,-0.041,-0.0344,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0;...
0,0,0,0,0,0,0,0,0,1,-1,0,0,0,0,-0.0804,-0.0389,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0;...
3,0,0,0,1,0,0,0,0,0,0,0,0,2.5,1.5,-8.6825,-9.2,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
0,0,0,2,1,0,1,0,0,0,0,-1,0,0,0,1.509,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0;...
0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,-0.0261,-0.024,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2,0,0;...
0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0;...
0,0,0,0,1,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.50;...
-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.452,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
0,0,-1,0,0,0,0,1,1,0,1,0,0,0,0,0.0082,-0.0479,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,-1,0;...
0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,-0.0873,-0.0449,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0;...
0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,-0.0807,-0.0719,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0;...
0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0;...
0,0,0,0,1,-1,-1,0,0,0,0,0,0,0,0,0.4445,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
2,-1,0,0,0,0,0,0,0,0,0,-1,-1,0,0,0.6391,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1;...
0,0,0,2,1,1,0,0,-1,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.50;...
0,0,0,0,0,0,0,1,-1,1,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0;...
0,0,0,0,0,0,0,0,0,0,0,0,0,-0.50,-0.50,-0.1002,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,-1,0,0;...
0,0,0,-1,0,1,0,0,0,0,1,0,0,0,0,-0.4270,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
2,-1,-1,-1,0,0,1,0,0,0,0,0,0,0,0,0.2085,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,-0.07130,-0.1000,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0;];
%     disp(leak)
%     disp(Aeq(11,15))
%     disp(Aeq(19,14))
    Aeq(11,15) = -leak;
    Aeq(19,14) = -leak;
    [met,rxns] = size(Aeq);
    beq = zeros(met,1);
    %% Reversible reactions
    rev = [2; 3; 8; 9; 21; 22; 24; 26; 27; 29; 33];

    %% Kinetic Reactions / Measured
    measindex = vmeas_index;
    vexp = vmeas;
%     save flux_measured measindex vexp;


    % Default UB, LB
    lb([1:rxns],1) = 0;
    ub([1:rxns],1) = inf;

    % Set strict for only reactions 10:13?
    strict = 10:13;
    lb(vmeas_index(strict)) = vmeas(strict) - abs(vmeas(strict))*err;
    ub(vmeas_index(strict)) = vmeas(strict) + abs(vmeas(strict))*err;

    % Set open LB for reversible
    lb(rev,1) = -inf;

    % Options
    oldoptions = optimset('fmincon');
    newoptions = optimset(oldoptions,'Display','off','LargeScale','off');

    % Initial values
    clear x0 x fval;
    x0 = rand(rxns,1)*100;

%% Test optimization structure
    % Making this into an anonymous function sped up excecution ~25x
    f = @(z) sum((   ((z(measindex) - vexp) ./ abs(vexp))  *100).^2);




    % Optimization command
    x = fmincon(f,x0,A,b,Aeq,beq,lb,ub,[],newoptions); %Jay- Why would they take this approach?
%     disp(x)
%     disp(vexp)

%     save workspace_flux;
%end
end