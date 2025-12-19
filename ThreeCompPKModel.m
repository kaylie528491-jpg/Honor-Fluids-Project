
classdef ThreeCompPKModel < handle
    % ThreeCompPKModel (IV bolus only, with binding; no explicit Jacobian)
    %
    % Three-compartment, amount-based PK model for Bob immediately post-flight.
    % Compartments: Plasma (central), Muscle (peripheral), Rest (peripheral).
    %
    % Dosing: SINGLE IV bolus at t = 0 (Dose mg injected into plasma).
    % Binding: Clearance and intercompartment exchange act on FREE drug only.
    %
    % Usage:
    %   PK = ThreeCompPKModel;         % construct with post-flight defaults
    %   PK.Dose = 100;                 % mg (optional change)
    %   PK.fu_p = 0.3;                 % fraction unbound in plasma (optional)
    %   PK.plot;                       % solve and plot Cp, Cm, Cr
    %   fprintf('AUC=%.2f  Cmax=%.2f  Vss=%.2f\n', PK.AUC, PK.Cmax, PK.Vss);

    %% === Changeable properties (post-flight Bob defaults) ===
    properties
        % Volumes L
        Vp  = 3.0 * (1 - 0.12);   % Plasma vol ~12% lower
        Vm  = 25.0 * (1 - 0.30);  % Muscle vol ~30% lower (atrophy)
        Vr  = 30.0;               % Rest-of-body vol

        % Intercompartmental clearances L/hr
        Qpm = 20 * 0.85;          % Plasma<->Muscle exchange (~15% reduction)
        Qpr = 15;                 % Plasma<->Rest exchange

        % Systemic clearance L/h (acts on free Cp)
        CL  = 5;

        % Partition coefficients (unitless)
        Kpm = 2.0;                % Muscle:plasma
        Kpr = 1.2;                % Rest:plasma

        % Binding (fractions unbound)
        fu_p = 0.4;               % plasma fraction unbound (edit per drug/physiology)
        fu_m = 1.0;               % muscle free fraction (often ~1 in interstitium)
        fu_r = 1.0;               % rest free fraction

        % Single bolus dose & simulation grid
        Dose    = 100;                       % mg (entire dose at t=0 into plasma)
        t_array = linspace(0,12,1200);       % hours
        fignum  = 10;                        % figure number
    end

    %% Dependent (computed) properties
    properties (Dependent)
        Numerical   % [N x 3] total concentrations [Cp, Cm, Cr]
        AUC         % mg*h/L (plasma total concentration)
        Cmax        % mg/L    (plasma total concentration)
        Vss         % L       (steady-state distribution volume; total-based)
    end

    %% Numerical solution (ode45)
    methods
        function Traj = get.Numerical(obj)
            % Amount-based ODE RHS
            fdot = @(t, A) obj.rhsAmounts(t, A);

            % Initial amounts: IV bolus in plasma at t=0
            A0 = [obj.Dose; 0; 0];   % [Ap(0), Am(0), Ar(0)] in mg

            % Let MATLAB estimate derivatives internally (no Jacobian provided)
            opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-9);
            [~, A] = ode45(fdot, obj.t_array, A0, opts);

            % Convert amounts to TOTAL concentrations [mg/L]
            Cp = A(:,1) / obj.Vp;
            Cm = A(:,2) / (obj.Vm * obj.Kpm);
            Cr = A(:,3) / (obj.Vr * obj.Kpr);

            Traj = [Cp, Cm, Cr];
        end

        function val = get.AUC(obj)
            C  = obj.Numerical(:,1);              % plasma total
            dt = obj.t_array(2) - obj.t_array(1);
            val = sum(C(1:end-1) + C(2:end)) * 0.5 * dt;
        end

        function val = get.Cmax(obj)
            C = obj.Numerical(:,1);
            val = max(C);
        end

        function v = get.Vss(obj)
            v = obj.Vp + obj.Kpm*obj.Vm + obj.Kpr*obj.Vr;
        end
    end

    %% ODE RHS (amounts), with binding
    methods
        function dA = rhsAmounts(obj, ~, A)
            % State vector: A = [Ap; Am; Ar] in mg
            Ap = A(1); Am = A(2); Ar = A(3);

            % TOTAL concentrations
            Cp_tot = Ap / obj.Vp;
            Cm_tot = Am / (obj.Vm * obj.Kpm);
            Cr_tot = Ar / (obj.Vr * obj.Kpr);

            % FREE concentrations (drive clearance & exchange)
            Cp_free = obj.fu_p * Cp_tot;
            Cm_free = obj.fu_m * Cm_tot;
            Cr_free = obj.fu_r * Cr_tot;

            % Intercompartment fluxes [mg/h] (FREE gradients)
            Jpm = obj.Qpm * (Cp_free - Cm_free);
            Jpr = obj.Qpr * (Cp_free - Cr_free);

            % Elimination on FREE plasma
            elim = obj.CL * Cp_free;

            % Amount ODEs [mg/h]
            dAp = -elim - Jpm - Jpr;
            dAm =  Jpm;
            dAr =  Jpr;

            dA  = [dAp; dAm; dAr];
        end
    end

    %% Plotting
    methods
        function plot(obj)
            C = obj.Numerical;   % TOTAL concentrations [Cp, Cm, Cr]
            Cp = C(:,1); Cm = C(:,2); Cr = C(:,3);

            % Optional: free plasma overlay
            Cp_free = obj.fu_p * Cp;

            figure(obj.fignum); clf;
            plot(obj.t_array, Cp,      'b-',  'LineWidth',2); hold on;
            plot(obj.t_array, Cm,      'r--', 'LineWidth',2);
            plot(obj.t_array, Cr,      'k-.', 'LineWidth',2);
            plot(obj.t_array, Cp_free, 'c:',  'LineWidth',1.5); % optional free plasma

            set(gca,'FontSize',14);
            xlabel('time (hr)'); ylabel('concentration (mg/L)');
            title('Bob (post-flight): Drug Concentrations in the body over time');
            legend({'Plasma (total)','Muscle (total)','Rest (total)','Plasma (free)'}, 'Location','best');
            grid on;
        end
    end
end