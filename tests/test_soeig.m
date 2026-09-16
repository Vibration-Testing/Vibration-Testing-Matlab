function test_soeig
%TEST_SOEIG Smoke test verifying SOEIG recovers known natural frequencies
% for a simple, analytically-known 2-DOF mass-spring system.
%
% Run with: test_soeig

% M = I, K = [2 -1; -1 2] has eigenvalues 1 and 3 (since M is identity,
% the eigenvalues of K itself are the squared natural frequencies in
% rad/s).
M = eye(2);
K = [2 -1; -1 2];

expected_fr = sort(sqrt([1 3])/(2*pi));
expected_fr = expected_fr(:);

[fr, ms] = soeig(M, K);
fr = sort(fr(:));

assert(numel(fr) == 2, 'soeig should return 2 natural frequencies for a 2-DOF system');
assert(max(abs(fr - expected_fr)) < 1e-10, ...
    'soeig natural frequencies do not match expected analytical values');
assert(isequal(size(ms), [2 2]), 'soeig should return a 2x2 mode shape matrix');

disp('test_soeig passed.');
end
