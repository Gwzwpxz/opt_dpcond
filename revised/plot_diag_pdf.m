clear; clc; close all;

l = load('reduce_left.mat'); l = l.reduce;
r = load('reduce_right.mat'); r = r.reduce;
lr = load('reduce_lr.mat'); lr = lr.reduce;

l = - log10(1 - l); l = l(l ~= inf);
r = - log10(1 - r);
lr = - log10(1 - lr);

% colormap(brewermap([], 'GnBu'));
% colormap(brewermap([],'RdYlGn'));
nhist(l, 'color', 'parula', 'binfactor', 3, 'LineWidth', 1);
yticklabels(round(yticks / length(l), 2));
ylabel("Percentage");
xlabel("Right log(cond before / cond after)");
% xlabel("Reduction in Right condition number");
set(gca, 'FontSize', 16);
tightfig;
saveas(gca, "right.pdf");

close all;

nhist(r, 'color', 'pink', 'binfactor', 3, 'LineWidth', 1);
yticklabels(round(yticks / length(r), 2));
ylabel("Percentage");
% xlabel("Reduction in Left condition number");
xlabel("Left log(cond before / cond after)");
set(gca, 'FontSize', 16);
tightfig;
saveas(gca, "left.pdf");

close all;

nhist(lr, 'color', 'lines', 'binfactor', 3, 'LineWidth', 1);
yticklabels(round(yticks / length(lr), 2));
ylabel("Percentage");
% xlabel("Reduction in Two-sided condition number");
xlabel("Two sided log(cond before / cond after)");
set(gca, 'FontSize', 16);
tightfig;
saveas(gca, "two.pdf");

close all;
subplot(3, 1, 1);
nhist(l, 'color', 'parula', 'binfactor', 3, 'LineWidth', 1, 'legend', "Right");
ylabel("");
set(gca, 'FontSize', 16);
yticklabels(round(yticks / length(l), 2));

subplot(3, 1, 2);
nhist(r, 'color', 'pink', 'binfactor', 3, 'LineWidth', 1, 'legend', "Left");
yticklabels(round(yticks / length(r), 2));
set(gca, 'FontSize', 16);
ylabel("Percentage");

subplot(3, 1, 3);
nhist(lr, 'color', 'lines', 'binfactor', 3, 'LineWidth', 1, 'legend', "Two-sided");
yticklabels(round(yticks / length(r), 2));
ylabel("");
xlabel("log(cond before / cond after)");
% xlabel("Reduction in condition number");

tightfig;
saveas(gca, "all.pdf");
