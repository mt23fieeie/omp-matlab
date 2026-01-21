function [x_hat, S, r_hist] = omp(A, b, K, tol, normalize_cols)
%OMP Orthogonal Matching Pursuit (Fig 3.1 template)
%
% Inputs:
%   A : (n x m) dictionary / measurement matrix
%   b : (n x 1) signal
%   K : max sparsity / max iterations
%   tol : stopping threshold on relative residual norm (optional)
%   normalize_cols : true/false, whether to normalize columns of A (optional)
%
% Outputs:
%   x_hat  : (m x 1) sparse coefficient estimate
%   S      : support set (selected indices)
%   r_hist : residual norm history ||r_k||_2
%
% Notes:
%   - Fig 3.1: each iteration selects one atom, then solves LS on selected atoms,
%              then updates residual.
%   - If normalize_cols = true, selection is based on normalized inner products
%     and final coefficients are mapped back to original A scaling.

    if nargin < 4 || isempty(tol) % nargin: parameter number
        tol = 0; % 0 means "do not stop early by tol"
    end
    if nargin < 5 || isempty(normalize_cols)
        normalize_cols = true;
    end

    [n, m] = size(A);
    assert(isvector(b) && length(b) == n, "b must be an (n x 1) vector."); % matched dimension
    b = b(:);% change to column vector

    % --- normalize columns  (Sec 3.1.4 Thm 3.1) ---
    col_norms = ones(1, m);
    A_work = A;
    if normalize_cols
        col_norms = vecnorm(A, 2, 1);          % 1 x m
        col_norms(col_norms == 0) = 1;         % avoid division by zero
        A_work = A ./ col_norms;               % normalize each column
    end

    % Fig 3.1: Initialize
    S = [];                                   % support set S_0 = empty
    r = b;                                    % residual r_0 = b
    r_hist = zeros(K, 1);

    % We will store the LS solution on support for each iteration
    xS = [];

    for k = 1:K
        % Fig 3.1: Select atom index j0
        projOntoOneVec = abs(A_work' * r);
        % Geometrically, project the r onto each column vector of A_work
        
        % m x 1, |<a_j, r>|
        projOntoOneVec(S) = 0;       
        % set those value with used index as zero               
        % do not re-select already chosen atoms

        [~, j0] = max(projOntoOneVec); % collect the column index

        % Fig 3.1: Update support set S_k = S_{k-1} U {j0}
        S = [S, j0];

        % Fig 3.1: Solve LS on selected atoms: min ||A_S x - b||_2
        AS = A_work(:, S);                    % n x k
        xS = AS \ b;             
        % (k x 1) least squares (projection), xS is the coordinate in the
        % projected subspace

        % Fig 3.1: Update residual r_k = b - A_S xS
        r = b - AS * xS;
        r_hist(k) = norm(r, 2);%l2 norm

        % Stopping Machanism: 
        % Each k, after applying Greedy algorithm, we check the overall
        % effect
        if tol > 0
            if r_hist(k) <= tol * norm(b,2)
                r_hist = r_hist(1:k);
                break;
            end
        end
    end

    % Build full x_hat (m x 1)
    x_hat = zeros(m, 1);
    if ~isempty(S)
        if normalize_cols
            % If we used A_work = A ./ col_norms, then:
            % b ≈ sum (A(:,j)/norm_j) * xS(j) = sum A(:,j) * (xS(j)/norm_j)
            % = sum A(:,j) * x_hat(S)
            x_hat(S) = xS ./ col_norms(S).'; 
        else
            x_hat(S) = xS;
        end
    end
end
