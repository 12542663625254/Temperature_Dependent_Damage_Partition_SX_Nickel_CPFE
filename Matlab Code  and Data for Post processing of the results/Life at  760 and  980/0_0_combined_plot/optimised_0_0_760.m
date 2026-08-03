clc; clear; close all;
% =========================================================================
%  STAGE 1 - CALIBRATE THE FATIGUE EXPONENT b ON PURE-FATIGUE DATA
%
%     N_cf = 1/(Dc^a + Df^b + M*(Dc*Df)^d)   with Dc = 0   ->   N = 1/Df^b
%
%  Only b is identifiable here (a, d, M, C, Qd multiply Dc = 0).
%  Carry the fitted b into Stage 2 (60/0 and 30/30 dwell data) to fit a,d,M.
% =========================================================================

%% 1. PURE-FATIGUE CASES (760 C = 1033 K)
cases(1).filename='cyclic_2_4.csv'; cases(1).epsilon_a=0.012; cases(1).exp_life=30;   cases(1).fw=[12 60];
cases(2).filename='cyclic_2_0.csv'; cases(2).epsilon_a=0.010; cases(2).exp_life=290;  cases(2).fw=[10 50];
cases(3).filename='cyclic_1_8.csv'; cases(3).epsilon_a=0.009; cases(3).exp_life=850;  cases(3).fw=[9 45];
cases(4).filename='cyclic_1_6.csv'; cases(4).epsilon_a=0.008; cases(4).exp_life=3557; cases(4).fw=[8 40];
exp_lives=[cases.exp_life];  N=numel(cases);

%% 2. FIXED Df-MODEL PARAMETERS (table values, same as the rest of the model)
Sg=0.38; Sc=0.9*Sg; Dfc=0.9; D0=0;
% (n1,Q,R,T,B1 are irrelevant here because there is no creep window -> Dc=0)

%% 3. COMPUTE Df FOR EACH CASE
Df_vals=zeros(1,N);
fprintf('Computing pure-fatigue Df...\n');
for k=1:N
    if ~isfile(cases(k).filename), error('File %s not found!',cases(k).filename); end
    raw=readmatrix(cases(k).filename);
    t=raw(:,1); s=raw(:,25); dSf=0;
    for i=2:numel(t)
        if (t(i)>cases(k).fw(1) && t(i)<=cases(k).fw(2)), dSf=dSf+s(i)*(t(i)-t(i-1)); end
    end
    if dSf>=Sg, Df=1.0;
    else, Df=D0 + ((Dfc-D0)/log(1-Sc/Sg))*log(max(1e-10,1-dSf/Sg)); end
    Df_vals(k)=max(1e-12,Df);
end

%% 4. FIT b   (N = 1/Df^b  ->  log10 N = -b*log10 Df, one parameter)
%   closed form (regression through origin) + bounded search for safety
x = log10(Df_vals);  y = log10(exp_lives);
b_closed = -sum(x.*y)/sum(x.^2);
obj = @(b) sqrt(mean((y - (-b.*x)).^2));
b = fminbnd(obj, 0.2, 4);
rmse = obj(b);

%% 5. PREDICT, REPORT, PLOT
pred = 1 ./ (Df_vals.^b);
err  = abs(pred-exp_lives)./exp_lives*100;

fprintf('\n=============== STAGE 1: FATIGUE EXPONENT b ===============\n');
fprintf(' b (fitted)        = %.4f\n', b);
fprintf(' b (closed-form)   = %.4f   (check)\n', b_closed);
fprintf(' log10-RMSE        = %.4f\n', rmse);
fprintf(' (a, d, M, C, Qd are NOT determined by pure-fatigue data)\n');
fprintf('==========================================================\n');
disp(table((1:N)',[cases.epsilon_a]',exp_lives',pred',err',Df_vals', ...
    'VariableNames',{'Case','Strain','Exp','Pred','Err_Pct','Df'}));

figure('Color','w','Position',[100 100 560 540]); hold on; box on;
allv=[exp_lives pred]; allv=allv(allv>0);
lim=[10^floor(log10(min(allv))) 10^ceil(log10(max(allv)))];
plot(lim,lim,'k-','LineWidth',1.4); plot(lim,2*lim,'k--'); plot(lim,0.5*lim,'k--');
scatter(exp_lives,pred,90,'filled','MarkerFaceColor',[0.10 0.45 0.80],'MarkerEdgeColor','k');
for k=1:N, text(exp_lives(k)*1.07,pred(k),sprintf('%.3f',cases(k).epsilon_a),'FontSize',9); end
set(gca,'XScale','log','YScale','log','FontSize',12); xlim(lim); ylim(lim); axis square; grid on;
xlabel('Experimental life (cycles)'); ylabel('Predicted life  N = 1/Df^b');
title(sprintf('Stage 1: pure fatigue, b = %.3f',b));
legend({'Perfect','x2 / /2',''},'Location','northwest');
saveas(gcf,'stage1_b_fit.png');
