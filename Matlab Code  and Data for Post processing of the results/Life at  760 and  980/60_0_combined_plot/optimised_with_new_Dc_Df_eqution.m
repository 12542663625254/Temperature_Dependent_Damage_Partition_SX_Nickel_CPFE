clc; clear; close all;
% =========================================================================
%  ROBUSTNESS TEST: apply 30/30-CALIBRATED interaction law to 60/0 DATA
%
%     N_cf = 1 / ( Dc^a + Df^b + M*(Dc*Df)^d )
%
%  Constants below are IMPORTED from the 30/30 fit and are NOT refitted.
%  Both datasets are at 760 C = 1033 K, so C and Qd do not enter here -
%  only a, b, d, M are used.  This is a pure generalisation check.
% =========================================================================

%% ---- CONSTANTS IMPORTED FROM THE 30/30 FIT (do NOT refit) ----
a = 0.9634;
b = 0.8807;
d = 2.9361;          % NOTE: hit the upper bound -> cubic, very scale-sensitive
M = 2.03811e8;       % lumped, tuned to the 30/30 Dc/Df SCALE

%% ---- B1 CHOICE (MUST match the recipe used during the 30/30 fit) ----
%   'table'     : B1 = 2.04e3 - 2.11e4*ea - 1.23e7*ea^2   (used in your fit)  <-- default
%   'sixtyzero' : B1 = 389501 - 630849*ea + 257797*ea^2   (your 60/0 script)
B1_MODE = 'table';
% WARNING: changing this changes the Dc scale by ~1e4. Because d=3 (cubic),
% an inconsistent B1 will make the transfer look wrong for the wrong reason.

%% 1. THE 60/0 CASES (760 C = 1033 K)
cases(1).filename='Fatigue_60_0_2_4.csv';     cases(1).epsilon_a=0.012; cases(1).exp_life=5;    cases(1).cw=[552 612]; cases(1).fw=[612 660];
cases(2).filename='Fatigue_60_0_0_8_perc.csv';cases(2).epsilon_a=0.008; cases(2).exp_life=600;  cases(2).cw=[100 160]; cases(2).fw=[160 200];
cases(3).filename='Fatigue_60_0_1_8.csv';     cases(3).epsilon_a=0.009; cases(3).exp_life=297;  cases(3).cw=[297 357]; cases(3).fw=[357 393];
cases(4).filename='Fatigue_60_0_2_0.csv';     cases(4).epsilon_a=0.010; cases(4).exp_life=100;  cases(4).cw=[310 370]; cases(4).fw=[370 410];
cases(5).filename='Fatigue_60_0_0_07.csv';    cases(5).epsilon_a=0.007; cases(5).exp_life=6600; cases(5).cw=[183 243]; cases(5).fw=[243 271];
exp_lives=[cases.exp_life];

%% 2. FIXED DAMAGE-MODEL PARAMETERS (same as the 30/30 model)
n1 = 0.6;  Q = 6.97e-19;  R = 8.314;  T = 1033;
Sg = 0.38; Sc = 0.9*Sg;   Dfc = 0.9;  D0 = 0;

%% 3. COMPUTE Dc, Df FOR EACH 60/0 CASE
Dc_vals=zeros(1,5); Df_vals=zeros(1,5);
fprintf('Computing 60/0 damages (B1_MODE = %s)...\n',B1_MODE);
for k=1:5
    ea = cases(k).epsilon_a;
    switch B1_MODE
        case 'table',     B1 = 2.04e3 - 2.11e4*ea - 1.23e7*ea^2;
        case 'sixtyzero', B1 = 389501 - 630849*ea + 257797*ea^2;
        otherwise, error('Unknown B1_MODE');
    end
    if ~isfile(cases(k).filename), error('File %s not found!',cases(k).filename); end
    raw = readmatrix(cases(k).filename);
    [Dc,dSf] = case_damage(raw(:,1), raw(:,25), B1, n1, Q, R, T, cases(k).cw, cases(k).fw);
    if dSf >= Sg, Df = 1.0;
    else, Df = D0 + ((Dfc-D0)/log(1-Sc/Sg))*log(max(1e-10,1-dSf/Sg)); end
    Dc_vals(k)=max(0,Dc); Df_vals(k)=max(0,Df);
end

%% 4. PREDICT (NO REFIT) - interaction law + linear rule for reference
pred_int = 1 ./ (Dc_vals.^a + Df_vals.^b + M*(Dc_vals.*Df_vals).^d);
pred_lin = 1 ./ (Dc_vals + Df_vals);
err_int  = abs(pred_int - exp_lives)./exp_lives*100;
err_lin  = abs(pred_lin - exp_lives)./exp_lives*100;
rmse_int = sqrt(mean((log10(exp_lives)-log10(max(1e-12,pred_int))).^2));
rmse_lin = sqrt(mean((log10(exp_lives)-log10(max(1e-12,pred_lin))).^2));
within2  = mean(pred_int>=exp_lives/2 & pred_int<=exp_lives*2)*100;
within3  = mean(pred_int>=exp_lives/3 & pred_int<=exp_lives*3)*100;

%% 5. REPORT
fprintf('\n========= ROBUSTNESS ON 60/0 (constants from 30/30 fit) =========\n');
fprintf(' a=%.4f  b=%.4f  d=%.4f  M=%.6g   (imported, not refitted)\n',a,b,d,M);
fprintf(' interaction-law  log10-RMSE = %.4f   (training value was 0.0128)\n',rmse_int);
fprintf(' linear-rule      log10-RMSE = %.4f\n',rmse_lin);
fprintf(' interaction predictions within x2: %.0f%%   within x3: %.0f%%\n',within2,within3);
fprintf('=================================================================\n');
disp(table((1:5)',[cases.epsilon_a]',exp_lives',pred_int',err_int',pred_lin',err_lin', ...
    Dc_vals',Df_vals', ...
    'VariableNames',{'Case','Strain','Exp','Pred_Int','ErrInt_Pct','Pred_Lin','ErrLin_Pct','Dc','Df'}));

%% 6. PLOT: PREDICTED vs EXPERIMENTAL (log-log, x2 & x3 bands)
figure('Color','w','Position',[100 100 560 540]); hold on; box on;
allv = [exp_lives pred_int pred_lin]; allv = allv(allv>0);
lim  = [10^floor(log10(min(allv))) 10^ceil(log10(max(allv)))];
plot(lim,lim,'k-','LineWidth',1.4);
plot(lim,2*lim,'k--'); plot(lim,0.5*lim,'k--');
plot(lim,3*lim,'b-.'); plot(lim,lim/3,'b-.');
hI=scatter(exp_lives,pred_int,90,'filled','MarkerFaceColor',[0.10 0.45 0.80],'MarkerEdgeColor','k');
hL=scatter(exp_lives,pred_lin,90,'d','MarkerEdgeColor','r','LineWidth',1.4);
for k=1:5, text(exp_lives(k)*1.07,pred_int(k),sprintf('%.3f',cases(k).epsilon_a),'FontSize',9); end
set(gca,'XScale','log','YScale','log','FontSize',12); xlim(lim); ylim(lim); axis square; grid on;
xlabel('Experimental life  (cycles)'); ylabel('Predicted life  (cycles)');
title(sprintf('30/30 constants on 60/0 data  (B1: %s)',B1_MODE));
legend([hI hL],{'Interaction (imported)','Linear rule'},'Location','northwest');
saveas(gcf,'robustness_60_0.png');

%% ===================== HELPERS ==========================================
function [Dc,dSf] = case_damage(t,s_dot_vec,B1,n1,Q,R,T,cw,fw)
    phi = B1*exp(-Q/(R*T)); Dc=0; dSf=0;
    for i=2:numel(t)
        dt=t(i)-t(i-1); s=s_dot_vec(i); tc=t(i);
        if (tc>cw(1)&&tc<=cw(2)) && s>0
            Dc = Dc + (1/phi)*(s)^(1-n1)*dt;
        elseif (tc>fw(1)&&tc<=fw(2))
            dSf = dSf + s*dt;
        end
    end
end