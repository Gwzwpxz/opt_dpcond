% Test the random sampling method
clear; clc; close all;

perturb = 0.0; matname = "imdb";
load imdb-opt.mat


rng(1);
X = data.X;
[m, n] = size(X);

nratio = 10;
ratios = logspace(-3, 0, nratio);
CovX = (data.M) / m;
nrmdiff = zeros(nratio, 1);
conds = zeros(nratio, 1);
Ddiag = diag(diag(data.M));
cond_diag = cond(sqrt(Ddiag) \ (sqrt(Ddiag) \ full(data.M))');

idx = 1;
for ratio = ratios
    msample = floor(m * ratio);
    randidx = randi(m, msample, 1);
    Xsample = X(randidx, :);
    CovXsample = (Xsample' * Xsample + speye(n) * perturb) / msample;
    nrmdiff(idx) = norm(CovXsample - CovX, 'fro');
    D = getcvxdiag(CovXsample, "L");
    conds(idx) = cond(sqrt(D) \ (sqrt(D) \ full(data.M))');
    idx = idx + 1;
end % End for

subplot(1, 2, 1);
semilogx(ratios, nrmdiff, "Marker", "*", "Color", "red", "MarkerSize", 10, "LineWidth", 3);
xlabel("m / M");
ylabel('|| X^T X - X_s ^T X_s ||');
set(gca, "FontSize", 12);
set(gca,'FontWeight','bold')

subplot(1, 2, 2);
semilogx(ratios, conds, "Marker", "*", "Color", "blue", "MarkerSize", 10, "LineWidth", 3);
ylabel('Condition number');
xlabel("m / M");
set(gca, "FontSize", 12);
set(gca,'FontWeight','bold')

tightfig;

saveas(gca, matname + ".fig");
saveas(gca, matname + ".pdf");

