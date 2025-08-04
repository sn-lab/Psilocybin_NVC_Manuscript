function [BOLD, Xfull, Ffull, Qfull, Vfull] = integrateBOLD(BOLD, X, Q, F, V, Z, dt, N, rho, alpha, V0, k1, k2, k3, Gamma, K, Tau)
    % Integrate the Balloon-Windkessel model
    % Adapted from neurolib/models/bold/timeIntegration.py
    % https://github.com/neurolib-dev/neurolib
    
    EPS = 1e-120;  % Small constant to prevent underflow
    Xfull = nan(size(Z));
    Qfull = nan(size(Z));
    Vfull = nan(size(Z));
    Ffull = nan(size(Z));

    for i = 1:size(Z, 2)  % Loop over all timesteps
        for j = 1:N  % Loop over all areas
            X(j) = X(j) + dt * (Z(j, i) - K(j) * X(j) - Gamma(j) * (F(j) - 1));
            Q(j) = Q(j) + dt / Tau(j) * (F(j) / rho * (1 - (1 - rho)^(1 / F(j))) - Q(j) * V(j)^(1 / alpha - 1));
            V(j) = V(j) + dt / Tau(j) * (F(j) - V(j)^(1 / alpha));
            F(j) = F(j) + dt * X(j);

            % Ensure F doesn't go below EPS to avoid errors
            F(j) = max(F(j), EPS);

            % Calculate the BOLD signal
            BOLD(j, i) = V0 * (k1 * (1 - Q(j)) + k2 * (1 - Q(j) / V(j)) + k3 * (1 - V(j)));

            Xfull(j,i) = X(j);
            Qfull(j,i) = Q(j);
            Vfull(j,i) = V(j);
            Ffull(j,i) = F(j);
        end
    end
end