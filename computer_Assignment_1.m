clear all
% Let's take inputs of rate expression parameters
fprintf('\n\n             [CHE331] COMPUTER ASSIGNMENT\n\n');
disp('The rate expression is: -rA = f(Ca) = k1Ca/(1 + k2Ca^2)');
disp('We are taking positive values of rate for simplicity. i.e. rA = k1Ca/(1 + k2Ca^2)');
k1 = input('\nEnter the value of rate constant k1: ');
k2 = input('Enter the value of rate constant k2: ');
Ca_i = input('Enter the initial value of Ca: ');
Ca_f = input('Enter the final value of Ca: ');

% k1 = 32;
% k2 = 3;
% Ca_i = 0;
% Ca_f = 20;
% Ca_d = 2.0;


% Create range of ca_i to ca_f with step size 0.01
Ca = Ca_i:0.001:Ca_f;

% Let's write reaction rate (-rA) as a function of Ca
rate = k1*Ca./(1+k2*Ca.^2);

% plot Ca vs rate
plot(Ca,rate,'k','LineWidth',1.5);
hold on
grid on
xlabel('Ca') 
ylabel('rate (rA)')

%(a)Let's Calculate derivative (of rate w.r.t. Ca) at a point Ca_d
drate = diff(rate)./diff(Ca);
Ca_d = input('\n(a): Enter the point at which you want to find derivative[d(rate)/d(Ca)]: ');
% derivative value
dtv_value = drate(find(Ca>=Ca_d,1));
fprintf('     The value of the derivative d(rate)/d(Ca) at ca=%.2f is %.6f.\n\n',Ca_d,dtv_value);

%(b)Let's find the area under the curve using trapezoidal rule of finding area
disp('(b): Enter the range in which area under the curve has to be calculated- ');
Ca1 = input('     Enter Initial value of the range: ');
Ca2 = input('     Enter Final value of the range: ');
Ca_range = Ca1:0.01:Ca2;
rate_range = k1*Ca_range./(1+k2*Ca_range.^2);
area = trapz(Ca_range,rate_range);
fprintf('     The area under the curve for given the range is: %.6f.\n\n',area);

%(c)Let's find the minimum and the maximum functional values
disp('(c): Enter the range in which maximum and minimum functional value has to be calculated- ');
Ca1 = input('     Enter Initial value of the range: ');
Ca2 = input('     Enter Final value of the range: ');
Ca_range = Ca1:0.01:Ca2;
rate_range = k1*Ca_range./(1+k2*Ca_range.^2);
[max_rate, ind_max] = max(rate_range);
[min_rate, ind_min] = min(rate_range);
% Let's find corresponding Ca values
max_Ca = Ca_range(ind_max);
min_Ca = Ca_range(ind_min);
% Let's Print the same and show both points on plot
fprintf('     The minimum functional value in the given range is: %.6f at Ca=%.2f.\n',min_rate,min_Ca);
fprintf('     i.e. (%.2f,%.6f)\n',min_Ca,min_rate);
scatter(min_Ca,min_rate, 50, 'r');
fprintf('     The maximum functional value in the given range is: %.6f at Ca=%.2f.\n',max_rate,max_Ca);
fprintf('     i.e. (%.2f,%.6f)\n\n',max_Ca,max_rate);
scatter(max_Ca,max_rate, 50, 'r');

%(d)Let's draw a straight line between two points (Ca1,rate1) and (Ca2,rate2)
disp('(d): Enter two points from which the straight line has to be drawn- ');
Ca1 = input('     Enter Ca value of first point: ');
Ca2 = input('     Enter Ca value of second point: ');
x = [Ca1,Ca2];
y = [rate(find(Ca>=Ca1,1)),rate(find(Ca>=Ca2,1))];
plot(x,y,'b','LineWidth',1.5);
% Let's find the slope and intercept of the line
p = polyfit(x, y, 1);
slope = p(1);
intercept = p(2);
fprintf('     The slope of the drawn line is: %.6f and the intercept is: %.6f.\n\n',slope,intercept);

% (e)Let's search for the point where the slope of the tangent is the same as the given slope of line drawn.
slope = input('(e): Enter a value of the slope: '); 
min_err_Ca = 0;
min_err_rate =0;
min_err = abs(drate(1)-slope);

for i = 1:length(Ca)-1
    if min_err > abs(drate(i)-slope)
        min_err = abs(drate(i)-slope);
        min_err_Ca = Ca(i);
        min_err_rate = rate(i);
    end
end

intercept = min_err_rate - slope*min_err_Ca;
len = input('     Enter length (projected on x-axis) of the tangent line: ');
x = [min_err_Ca-len/2,min_err_Ca+len/2];
y = [(slope*(min_err_Ca-len/2)+intercept),(slope*(min_err_Ca+len/2)+intercept)];
plot(x,y,'g','LineWidth',1.5);
fprintf('     At the point (%.2f,%.6f) the slope of the tangent is the same as the given slope.\n',min_err_Ca,min_err_rate);
