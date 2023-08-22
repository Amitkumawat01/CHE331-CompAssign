Ca0 = [2 5 6 6 11 14 16 24];
Ca = [0.5 3 1 2 6 10 8 4];
tau = [30 1 50 8 4 20 20 4];

neg_rate_inv = tau./[Ca0 - Ca];
disp(neg_rate_inv);

% Fit cubic spline
pp = spline(Ca, neg_rate_inv);

Ca_interp = linspace(min(Ca), max(Ca), 1000);
neg_rate_inv_interp = ppval(pp, Ca_interp);

% Given data
v0 = 0.1;
Ca1 = 1;
Ca2 = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) A Single PFR

subplot(2,3,1)
% Plot original data points and spline
plot(Ca, neg_rate_inv, 'o', 'MarkerFaceColor','none', 'MarkerEdgeColor', 'b', 'MarkerSize', 8);
hold on;
plot(Ca_interp, neg_rate_inv_interp, 'r-', 'LineWidth', 2);

% labels and title
xlabel('Ca');
ylabel('-1/r_A');
title('A Single PFR');
grid on;
% Calculate area under the curve
area = integral(@(x) ppval(pp, x), Ca1, Ca2);

fprintf('(A) A Single PFR\n');
fprintf('Minimum volume required: %.4f\n', v0*area);

% Shade the area under the curve
x_fill = linspace(Ca1, Ca2, 100);
y_fill = ppval(pp, x_fill);
fill([Ca1, x_fill, Ca2], [0, y_fill, 0], 'g', 'FaceAlpha', 0.2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (b) A single CSTR

subplot(2,3,2)
% Plot original data points and spline
plot(Ca, neg_rate_inv, 'o', 'MarkerFaceColor','none', 'MarkerEdgeColor', 'b', 'MarkerSize', 8);
hold on;
plot(Ca_interp, neg_rate_inv_interp, 'r-', 'LineWidth', 2);

% labels and title
xlabel('Ca');
ylabel('-1/r_A');
title('A single CSTR');
grid on;

area = ppval(pp, Ca1)*(Ca2 - Ca1);
fprintf('\n(B) A single CSTR\n');
fprintf('Minimum volume required: %.4f\n', v0*area);

% Shade the area under the curve
x_fill = linspace(Ca1, Ca2, 100);
y_fill = ones(1, length(x_fill)) * ppval(pp, Ca1);
fill([Ca1, x_fill, Ca2], [0, y_fill, 0], 'g', 'FaceAlpha', 0.2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Two MFR 

subplot(2,3,3)
% Plot original data points and spline
plot(Ca, neg_rate_inv, 'o', 'MarkerFaceColor','none', 'MarkerEdgeColor', 'b', 'MarkerSize', 8);
hold on;
plot(Ca_interp, neg_rate_inv_interp, 'r-', 'LineWidth', 2);

% labels and title
xlabel('Ca');
ylabel('-1/r_A');
title('Two MFR');
grid on;

area = 100000;
Caf  = 0;

for Cai = Ca1:0.01:Ca2
    area1 = (Cai - Ca1)*ppval(pp, Ca1);
    area2 = (Ca2 - Cai)*ppval(pp, Cai);
    
    if area1 + area2 < area
        area = area1 + area2;
        Caf = Cai;
    end
end

fprintf('\n(C) Two MFR\n');
fprintf('Minimum volume required: %.4f\n', v0*area);

% Shade the area under the curve
x_fill = linspace(Ca1, Caf, 100);
y_fill = ones(1, length(x_fill)) * ppval(pp, Ca1);
fill([Ca1, x_fill, Caf], [0, y_fill, 0], 'g', 'FaceAlpha', 0.2);

% Shade the area under the curve
x_fill = linspace(Caf, Ca2, 100);
y_fill = ones(1, length(x_fill)) * ppval(pp, Caf);
fill([Caf, x_fill, Ca2], [0, y_fill, 0], 'g', 'FaceAlpha', 0.2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (D) PFR with MFR

subplot(2,3,4)
% Plot original data points and spline
plot(Ca, neg_rate_inv, 'o', 'MarkerFaceColor','none', 'MarkerEdgeColor', 'b', 'MarkerSize', 8);
hold on;
plot(Ca_interp, neg_rate_inv_interp, 'r-', 'LineWidth', 2);

% labels and title
xlabel('Ca');
ylabel('-1/r_A');
title('PFR with MFR');
grid on;

% Finding Ca at minima 
objective = @(x) ppval(pp, x);
initial_guess = 0.5;
minima_ca = fminbnd(objective, min(Ca), max(Ca));

area1 = integral(@(x) ppval(pp, x), Ca1, minima_ca);
area2 = ppval(pp, minima_ca)*(Ca2 - minima_ca);
area = area1 + area2;

fprintf('\n(D) PFR with MFR\n');
fprintf('Minimum volume required: %.4f\n', v0*area);

% Shade the area under the curve
x_fill = linspace(Ca1, minima_ca, 100);
y_fill = ppval(pp, x_fill);
fill([Ca1, x_fill, minima_ca], [0, y_fill, 0], 'g', 'FaceAlpha', 0.2);

% Shade the area under the curve
x_fill = linspace(minima_ca, Ca2, 100);
y_fill = ones(1, length(x_fill)) * ppval(pp, minima_ca);
fill([minima_ca, x_fill, Ca2], [0, y_fill, 0], 'b', 'FaceAlpha', 0.2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (E) A PFR with Recycle 
subplot(2,3,5)
% Plot original data points and spline
plot(Ca, neg_rate_inv, 'o', 'MarkerFaceColor','none', 'MarkerEdgeColor', 'b', 'MarkerSize', 8);
hold on;
plot(Ca_interp, neg_rate_inv_interp, 'r-', 'LineWidth', 2);

%labels and title
xlabel('Ca');
ylabel('-1/r_A');
title('A PFR with Recycle ');
grid on;

min_diff = 100000;
Caf  = 0;

for Cai = Ca1:0.01:Ca2
    LHS = ppval(pp,Cai);
    RHS = (integral(@(x) ppval(pp, x), Ca1, Cai))/(Cai-Ca1);
    
    if abs(LHS - RHS) < min_diff
        min_diff = abs(LHS - RHS);
        Caf = Cai;
    end
end
%disp(min_diff);

area = ppval(pp, Caf)*(Ca2-Ca1);
disp(Caf);
disp(ppval(pp, Caf));
fprintf('\n(E) A PFR with Recycle\n');
fprintf('Minimum volume required: %.4f\n', v0*area);

% Shade the area under the curve
x_fill = linspace(Ca1, Ca2, 100);
y_fill = ones(1, length(x_fill)) * ppval(pp, Caf);
fill([Ca1, x_fill, Ca2], [0, y_fill, 0], 'g', 'FaceAlpha', 0.2);