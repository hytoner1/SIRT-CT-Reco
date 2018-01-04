function x = simple_sirt(A, b)
% Input: sparse system matrix A, data b.
% Output: SIRT reconstruction x.
% x = zeros(d * d, 1);
[rows, cols] = size(A);
x = zeros(cols, 1);
C = sparse(1 : cols, 1 : cols, 1 ./ sum(A));
R = sparse(1 : rows, 1 : rows, 1 ./ sum(A'));
CATR = C * A' * R;
for i = 1 : 100
  x = x + CATR * (b - A * x);
  imagesc(reshape(x,100,100))
end

end