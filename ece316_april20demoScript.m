% Saving a simulation from Simulink

clear
clc

model_output = sim('ece316_april20demo.slx','ReturnWorkspaceOutputs','on');

y = model_output.yout{1}.Values.Data;
t = model_output.yout{1}.Values.Time;

figure(1)
plot(t,y)
set(gca,'Fontsize',16)
xlabel('time','Fontsize',16)
ylabel('H(s)','Fontsize',16)

%modeloutput = yout{1}.Values.Data;
%plot(modeloutput)

