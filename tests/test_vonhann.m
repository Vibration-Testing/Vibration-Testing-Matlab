function test_vonhann
%TEST_VONHANN Smoke test verifying VONHANN (and HANNWIN) produce a
% correctly shaped, symmetric window with values in [0, 1].
%
% Run with: test_vonhann

n = 64;
w = hannwin(n);

assert(numel(w) == n, 'hannwin should return a window of length n');
assert(all(w >= -1e-10) && all(w <= sqrt(8/3) + 1e-10), ...
    'hannwin values should lie within [0, sqrt(8/3)] (its energy-normalized range)');
assert(max(abs(w - flipud(w))) < 1e-10, 'hannwin window should be symmetric');
assert(abs(w(1)) < 0.05, 'hannwin should taper to (near) zero at the edges');

disp('test_vonhann passed.');
end
