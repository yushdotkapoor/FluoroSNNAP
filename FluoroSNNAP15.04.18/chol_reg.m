
function [R, flag] = chol_reg(A, triangle)

    eig_EE=eig(A);
    if all(eig_EE > 0)
        [R, flag] = chol(A, triangle); % inverse covariance square root
    else
        disp('Matrix is not positive definite, attempting to regularize');
        lambda = 0.01; % regularization parameter
        A_reg = A + lambda*eye(size(A));
        [R, flag] = chol(A_reg, triangle); % inverse covariance square root
    end