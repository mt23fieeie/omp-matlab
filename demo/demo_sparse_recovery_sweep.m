% demo_sparse_recovery_sweep.m
% ------------------------------------------------------------
% OMP demo (two versions in one script):
%
% Version A (single run):
%   - Run OMP once for a chosen noise level sigma.
%   - Save one residual plot and one .mat file.
%
% Version B (noise sweep):
%   - Run OMP for multiple sigma values with a fair comparison:
%     A and x0 are identical across sigmas (same rng seed),
%     only the additive noise changes with sigma.
%   - Print a summary table and save one residual plot per sigma.
%
% Notes:
%   - This script assumes src/omp.m is on MATLAB path.
%   - Outputs are saved into ./results/ (created if needed).
% ------------------------------------------------------------

clear; clc;

%% ============== Version A: Single run (pick one sigma) ==============
do_single_run = false;     % set false if you only want the sweep
sigma_single = 10;        % try: 0.01, 0.1, 1, 10
tol_single   = 1e-8;      % strict relative residual threshold
normalize_cols = true;    % normalize dictionary columns internally

if do_single_run
    rng(0); % Fix random seed for reproducibility

    % Problem dimensions
    n = 64; m = 128; K = 6;

    % Dictionary
    A = randn(n, m);

    % K-sparse ground truth
    x0 = zeros(m, 1);
    S0 = randperm(m, K);
    x0(S0) = randn(K, 1);

    % Noisy observation
    b = A * x0 + sigma_single * randn(n, 1);

    % Run OMP
    [x_hat, S_hat, r_hist] = omp(A, b, K, tol_single, normalize_cols);

    % Metrics
    hit = numel(intersect(S0, S_hat));
    rel_err_x = norm(x_hat - x0) / norm(x0);
    rel_resid = norm(b - A * x_hat) / norm(b);

    fprintf("=== Single run (sigma = %g) ===\n", sigma_single);
    fprintf("hit: %d/%d\n", hit, K);
    fprintf("rel err x: %.3e\n", rel_err_x);
    fprintf("rel residual: %.3e\n\n", rel_resid);

    % Save outputs
    if ~exist("results", "dir"), mkdir("results"); end
    figure; plot(r_hist, '-o'); grid on;
    xlabel('k'); ylabel('||r_k||_2');
    title(sprintf('OMP residual (sigma=%g)', sigma_single));
    exportgraphics(gcf, sprintf("results/residual_sigma%g.png", sigma_single));
    save(sprintf("results/run_sigma%g.mat", sigma_single), ...
        "n","m","K","sigma_single","S0","S_hat","x0","x_hat","r_hist","hit","rel_err_x","rel_resid");
    close;
end

%% ============== Version B: Sweep over sigmas (fair comparison) ==============
do_sweep = true;          % set false if you only want single run
sigmas = [0.01, 0.1, 1, 10];
tol_sweep = 0;            % use fixed K iterations in sweep (no early stop)

if do_sweep
    fprintf("\n%-8s %-8s %-12s %-12s\n","sigma","hit","rel_err_x","rel_resid");
    fprintf("----------------------------------------------------\n");

    for s = sigmas
        rng(0); % Ensure A and x0 are identical across sigmas (fair comparison)

        % Problem dimensions
        n = 64; m = 128; K = 6;

        % Dictionary
        A = randn(n, m);

        % K-sparse ground truth
        x0 = zeros(m, 1);
        S0 = randperm(m, K);
        x0(S0) = randn(K, 1);

        % Noisy observation
        b = A * x0 + s * randn(n, 1);

        % Run OMP
        [x_hat, S_hat, r_hist] = omp(A, b, K, tol_sweep, normalize_cols);

        % Metrics
        hit = numel(intersect(S0, S_hat));
        rel_err_x = norm(x_hat - x0) / norm(x0);
        rel_resid = norm(b - A * x_hat) / norm(b);

        fprintf("%-8g %-8d %-12.3e %-12.3e\n", s, hit, rel_err_x, rel_resid);

        % Save residual plot + mat log
        if ~exist("results", "dir"), mkdir("results"); end
        figure; plot(r_hist,'-o'); grid on;
        xlabel('k'); ylabel('||r_k||_2');
        title(sprintf('OMP residual (sigma=%g)', s));
        exportgraphics(gcf, sprintf("results/residual_sigma%g.png", s));
        save(sprintf("results/run_sigma%g.mat", s), ...
            "n","m","K","s","S0","S_hat","x0","x_hat","r_hist","hit","rel_err_x","rel_resid");
        close;
    end
end
