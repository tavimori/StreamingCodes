function result = simulate_all(alpha,beta,sim_length,channel_type) 
eps = linspace(0.01,0.2,2);  %Select the range of epsilon
l = length(eps);

no_rates = 6;
result = -1*ones(l,no_rates+5);

if channel_type == "GE"
    start_time = strcat('GE_',datestr(now,'mm-dd_HH-MM'));
elseif channel_type == "Fritchman"
    start_time = strcat('Fritchman_',datestr(now,'mm-dd_HH-MM'));
else
    fprintf("unknown channel model!\n");
end
%start_time = stract('GE_',datestr(now,'mm-dd_HH-MM'));

mkdir(start_time);
Filename = strcat(start_time,'/result.txt');
tic

% parfor i = 1:l %Replace by for if the system does not have multiple cores
for i = 1:l %Replace by for if the system does not have multiple cores
    fprintf('\n Running Simulation for epsilon = %f\n',eps(i))
    result(i,:) =  simulate_once(alpha,beta,sim_length,eps(i));
    dlmwrite(Filename,result(i,:),'delimiter','\t','-append');
    
end
toc

% unsorted_result = dlmread('fritchman_compiled.txt'); %%REPLACE IT BY THE NAME OF OUTPUT FILE
unsorted_result = dlmread(Filename); %%REPLACE IT BY THE NAME OF OUTPUT FILE
[~,idx] = sort(unsorted_result(:,1)); % sort just the first column
result = unsorted_result(idx,:);
figure('units','normalized','outerposition',[0 0 1 1])
%hold on
eps = result(:,1);
mds = result(:,10);
const_c = result(:,9);
mt = result(:,8);
fo_kh = result(:,7);
const_a = result(:,6);
uncoded = result(:,5);
th_uncoded = result(:,11);
%semilogy(eps,result(:,5),'-*r',eps,result(:,11),'-.c',eps,result(:,6),'-sg',eps,result(:,7),'-ob',eps,result(:,8),'-.m',eps,result(:,9),'-sc',eps,result(:,10),'--dk')
semilogy(eps,mds,eps,mt,eps,fo_kh,eps,const_a,eps,const_c); 
grid on
legend({'MDS','Martinian-Trott','Fong et al.','Construction A','Construction C'},'Orientation','horizontal','Location','northoutside','FontSize',16)
%legend({'Uncoded','Theoretical Uncoded','Construction A','Fong Khishti','Martinian-Trott','Construction C', 'MDS'},'Orientation','horizontal','Location','northoutside','FontSize',16)
xlabel('Epsilon')
ylabel('Packet Loss Rate')
% set(gcf,'Units','Inches');
% pos = get(gcf,'Position');
% set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

%MDS MT FONG A C 






















end