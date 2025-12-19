%% MuscleFVM_Bob  (Fluid-Control-Volume Balance Form)
%
% Finite-Volume (1D) solute transport in interstitial FLUID,
% driven by plasma concentration Cp(t) at x=0 via an exchange flux.
%
% Fluid-CV statement (per unit cross-sectional area):
%   For each fluid slab i of thickness Dx:
%     Accumulation in fluid slab = (solute flux in) - (solute flux out) - (uptake from fluid)
%
%   d/dt ( C_i * Dx ) = J_{i-1/2} - J_{i+1/2} - k_uptake * C_i * Dx
%   Divide by Dx:
%     dC_i/dt = (J_{i-1/2} - J_{i+1/2})/Dx - k_uptake * C_i
%
% Here the solute flux in the fluid is diffusive (Fick's law):
%   J = -D * dC/dx
%
% Boundary fluxes (fluid-side):
%   x=0 exchange (Robin-type):  J_{1/2} = k_exchange * ( Cp(t) - C(0,t) )
%   x=L no-flux:               J_{N+1/2} = 0
%
% So this is a FLUID BALANCE on concentration in each interstitial fluid element.

classdef MuscleFVM_Bob
    %% ===== Changeable properties =====
    properties
        L (1,1) {mustBePositive} = 2.0
        N (1,1) {mustBeInteger,mustBePositive} = 100
        D (1,1) {mustBePositive} = 3.6e-3
        k_uptake (1,1) {mustBeNonnegative} = 0.10
        k_exchange (1,1) {mustBeNonnegative} = 0.30
        t_array (1,:) double = linspace(0,12,1200)
        useImplicit (1,1) logical = false
    end

    %% ===== Derived (private helpers) =====
    properties (Access=private)
        Dx
        dt
        Nt
        A_mat
    end

    %% ===== Public API =====
    methods
        function obj = MuscleFVM_Bob(varargin)
            if ~isempty(varargin)
                p = inputParser;
                addParameter(p,'L',obj.L);
                addParameter(p,'N',obj.N);
                addParameter(p,'D',obj.D);
                addParameter(p,'k_uptake',obj.k_uptake);
                addParameter(p,'k_exchange',obj.k_exchange);
                addParameter(p,'t_array',obj.t_array);
                addParameter(p,'useImplicit',obj.useImplicit);
                parse(p,varargin{:});
                f = fieldnames(p.Results);
                for k=1:numel(f), obj.(f{k}) = p.Results.(f{k}); end
            end
            obj.Dx = obj.L / obj.N;
            obj.Nt = numel(obj.t_array);
            obj.dt = obj.t_array(2) - obj.t_array(1);
        end

        function Ct = Profiles(obj, Cp_time)
            % Return [Nt x N] interstitial-fluid concentration profile [mg/L].
            validateattributes(Cp_time, {'double'},{'vector'});
            if numel(Cp_time) ~= obj.Nt
                error('Cp_time length (%d) must equal t_array length (%d).', numel(Cp_time), obj.Nt);
            end

            % ---- Quick stability (explicit diffusion) check ----
            if ~obj.useImplicit
                cfl = obj.Dx^2/(2*obj.D);
                if obj.dt > cfl
                    warning('Explicit step may be unstable: dt=%.4g h > Dx^2/(2D)=%.4g h', obj.dt, cfl);
                end
            else
                obj = obj.buildImplicitSystem();
            end

            % Concentration in each interstitial-fluid CV cell
            Ct = zeros(obj.Nt, obj.N);

            for n = 1:obj.Nt-1
                C  = Ct(n,:).';     % current fluid concentration profile [N x 1]
                Cp = Cp_time(n);    % plasma concentration at current time

                if ~obj.useImplicit
                    % =========================================================
                    % EXPLICIT FLUID CONTROL-VOLUME BALANCE
                    %
                    % Define diffusive solute flux at internal faces:
                    %   J_{i+1/2} = -D (C_{i+1} - C_i)/Dx
                    % =========================================================
                    J = -obj.D * diff(C) / obj.Dx;   % size N-1, internal faces

                    % Fluid balance per cell:
                    %   dC_i/dt = (J_{i-1/2} - J_{i+1/2})/Dx - k_uptake*C_i
                    dC = zeros(obj.N,1);

                    % Interior cells (2..N-1)
                    dC(2:obj.N-1) = (J(1:obj.N-2) - J(2:obj.N-1))/obj.Dx ...
                                    - obj.k_uptake * C(2:obj.N-1);

                    % Left boundary face (x=0): exchange flux INTO fluid cell 1
                    %   J_{1/2} = k_exchange*(Cp - C1)
                    J_left = obj.k_exchange * (Cp - C(1));
                    dC(1)  = (J_left - J(1))/obj.Dx - obj.k_uptake * C(1);

                    % Right boundary face (x=L): no-flux (Neumann)
                    %   J_{N+1/2} = 0
                    dC(obj.N) = (J(obj.N-1) - 0.0)/obj.Dx - obj.k_uptake * C(obj.N);

                    % Time update (fluid concentration in each CV)
                    C_next = C + obj.dt * dC;
                    C_next = max(C_next, 0);
                    Ct(n+1,:) = C_next.';

                else
                    % ---------- Implicit Backward Euler ----------
                    rhs = C;
                    Cp_next = Cp_time(n+1);
                    rhs(1) = rhs(1) + obj.dt * (obj.k_exchange/obj.Dx) * Cp_next;

                    C_next = obj.A_mat \ rhs;
                    C_next = max(C_next, 0);
                    Ct(n+1,:) = C_next.';
                end
            end
        end

        function plotProfiles(obj, Cp_time, snap_hours)
            Ct = obj.Profiles(Cp_time);
            x  = linspace(0, obj.L, obj.N);

            figure; hold on;
            for k = 1:numel(snap_hours)
                [~, idx] = min(abs(obj.t_array - snap_hours(k)));
                plot(x, Ct(idx,:), 'LineWidth', 1.8, ...
                    'DisplayName', sprintf('t = %.1f h', obj.t_array(idx)));
            end
        end
    end
end
