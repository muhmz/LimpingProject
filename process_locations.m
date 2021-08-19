
t1 = [];
t2 = [];
t_all = [];
for i = 1:size(results.locs,1)
    
    
    gait_steps =  results.locs{i,1};
    
    difference = diff(gait_steps);
    
    tau1 = difference(1:2:length(difference));
    tau2 = difference(2:2:length(difference));
    
    
    t1 = [t1, tau1];
    t2 = [t2, tau2];
    
    t_all = [t_all, difference];
end 
%figure;

% Remove outliers
t1 = rmoutliers(t1,'percentiles',[10 90]);
t2 = rmoutliers(t2,'percentiles',[10 90]);
%t_all = rmoutliers(t_all,'percentiles',[10 90]);


figure;
j = [0:0.01:2];
[f1,xi1] = ksdensity(t1,j ,'Bandwidth',0.03);
plot(xi1,f1);

hold on
[f2,xi2] = ksdensity(t2,j, 'Bandwidth',0.03);
plot(xi2,f2);

% [f3,xi3] = ksdensity(t_all,j, 'Bandwidth',0.05);
% plot(xi3,f3);

grid on;
ax.FontSize = 10;
ylabel('Density (f)','Interpreter','latex', 'FontSize',16);
xlabel('X$_i$', 'Interpreter','latex', 'FontSize',16);

hold off

figure;
% scatter plot
g = [t1, t2];
l = [1:1:size(g,2)];
scatter(l,g)
xlim([0 length(l)]);
ylim([0 2]);