%addpath('.')
table = readtable('test.txt');
y = table.Glucose;
t = table.DateTime;
smoother_result = SmoothSMBGData(t,y,'outlierRemoval',1,'dynamicModel',2);
%Put smoother result into table as new column
table.SmoothedGlucose = smoother_result.y_smoothed_at_tout;
%To store table with smoothed glucose back to file, uncomment below
%writetable(table,'test_smoothed.txt')

%plots
plot(t,y,'r.','MarkerSize',20)
hold on
ol = smoother_result.outliers==1;
plot(t(ol),y(ol),'kx','MarkerSize',10)
plot(smoother_result.t_i,smoother_result.y_smoothed,'r-');
hold off
legend('Input glucose measurements','Outliers','Smoothed glucose','location','South')




      