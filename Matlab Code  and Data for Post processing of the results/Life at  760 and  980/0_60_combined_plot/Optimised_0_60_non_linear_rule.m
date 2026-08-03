clc; clear; close all;
% =========================================================================
%  ROBUSTNESS TEST on 0/60 (compressive-dwell) data, constants from 30/30
%
%     N_cf = 1/( Dc^a + Df^b + M*(Dc*Df)^d )
%
%  Same B1 formula and same T = 1033 K as the 30/30 calibration, so the
%  Dc scale is consistent -> a fair transfer test. C, Qd do not enter at
%  one temperature; only a, b, d, M are used. NO refitting.
% =========================================================================

%% ---- CONSTANTS IMPORTED FROM THE 30/30 STAGE-2 FIT (do NOT refit) ----
a = 0.9634;          % <-- PASTE a from the Stage-2 (30/30) run
b = 0.8807;       % fixed (Stage 1, pure fatigue)
d = 2.9361;          % <-- PASTE d from the Stage-2 (30/30) run
M = 1.27131e8;    % lumped, from the Stage-2 (30/30) run

%% 1. 0/60 CASES (760 C = 1033 K).  cw = creep window, fw = fatigue windows
cases(1).filename='Pin_2_4_60C.csv';     cases(1).epsilon_a=0.012; cases(1).exp_life=16;  cases(1).cw=[252 312]; cases(1).fw=[228 252;312 336];
cases(2).filename='Pin_2_0_60C.csv';     cases(2).epsilon_a=0.010; cases(2).exp_life=176; cases(2).cw=[330 390]; cases(2).fw=[310 330;390 410];
cases(3).filename='Pin_1_8_60C.csv';     cases(3).epsilon_a=0.009; cases(3).exp_life=318; cases(3).cw=[315 375]; cases(3).fw=[297 315;375 393];
cases(4).filename='Pin_1_6_60_Comp.csv'; cases(4).epsilon_a=0.008; cases(4).exp_life=534; cases(4).cw=[116 176]; cases(4).fw=[100 116;176 192];
exp_lives=[cases.exp_life]; N=numel(cases);

%% 2. FIXED DAMAGE-MODEL PARAMETERS (same as the 30/30 calibration)
n1=0.6; Q=6.97e-19; R=8.314; T=1033;
Sg=0.38; Sc=0.9*Sg; Dfc=0.9; D0=0;

if isnan(a)||isnan(d), error('Set a and d (from the Stage-2 30/30 fit) before running.'); end

%% 3. COMPUTE Dc, Df FOR EACH 0/60 CASE  (table B1, consistent with 30/30)
Dc_vals=zeros(1,N); Df_vals=zeros(1,N);
fprintf('Computing 0/60 damages...\n');
for k=1:N
    ea=cases(k).epsilon_a;
    B1 = 2.04e3 - 2.11e4*ea - 1.23e7*ea^2;        % table formula (matches calibration)
    if ~isfile(cases(k).filename), error('File %s not found!',cases(k).filename); end
    raw=readmatrix(cases(k).filename);
    [Dc,dSf]=damage(raw(:,1),raw(:,25),B1,n1,Q,R,T,cases(k).cw,cases(k).fw);
    if dSf>=Sg, Df=1.0;
    else, Df=D0+((Dfc-D0)/log(1-Sc/Sg))*log(max(1e-10,1-dSf/Sg)); end
    Dc_vals(k)=max(0,Dc); Df_vals(k)=max(0,Df);
end

%% 4. PREDICT (NO REFIT) - interaction + linear for reference
% --- individual denominator terms of the interaction law ---
term1_Dc  = Dc_vals.^a;                  % first term:  Dc^a
term2_Df  = Df_vals.^b;                  % second term: Df^b
term3_int = M*(Dc_vals.*Df_vals).^d;     % coupling term: M*(Dc*Df)^d
denom_int = term1_Dc + term2_Df + term3_int;

pred_int = 1 ./ denom_int;
pred_lin = 1 ./ (Dc_vals + Df_vals);
err_int  = abs(pred_int-exp_lives)./exp_lives*100;
err_lin  = abs(pred_lin-exp_lives)./exp_lives*100;
rmse_int = sqrt(mean((log10(exp_lives)-log10(max(1e-12,pred_int))).^2));
rmse_lin = sqrt(mean((log10(exp_lives)-log10(max(1e-12,pred_lin))).^2));
within2  = mean(pred_int>=exp_lives/2 & pred_int<=exp_lives*2)*100;
within3  = mean(pred_int>=exp_lives/3 & pred_int<=exp_lives*3)*100;

% --- fractional contribution of each term to the denominator ---
frac1 = term1_Dc  ./ denom_int * 100;
frac2 = term2_Df  ./ denom_int * 100;
frac3 = term3_int ./ denom_int * 100;

%% 5. REPORT
fprintf('\n======= ROBUSTNESS ON 0/60 (constants from 30/30 fit) =======\n');
fprintf(' a=%.4f  b=%.4f  d=%.4f  M=%.6g   (imported, not refitted)\n',a,b,d,M);
fprintf(' interaction-law  log10-RMSE = %.4f\n',rmse_int);
fprintf(' linear-rule      log10-RMSE = %.4f\n',rmse_lin);
fprintf(' interaction within x2: %.0f%%   within x3: %.0f%%\n',within2,within3);
fprintf('=============================================================\n');
disp(table((1:N)',[cases.epsilon_a]',exp_lives',pred_int',err_int',pred_lin',err_lin',Dc_vals',Df_vals', ...
    'VariableNames',{'Case','Strain','Exp','Pred_Int','ErrInt_Pct','Pred_Lin','ErrLin_Pct','Dc','Df'}));

% --- NEW: breakdown of the three interaction-law terms ---
fprintf('\n---------- INTERACTION-LAW TERM BREAKDOWN  (denominator of N_cf) ----------\n');
fprintf('  N_cf = 1 / ( Dc^a  +  Df^b  +  M*(Dc*Df)^d )\n\n');
disp(table((1:N)',[cases.epsilon_a]', ...
    term1_Dc', term2_Df', term3_int', denom_int', ...
    frac1', frac2', frac3', ...
    'VariableNames',{'Case','Strain', ...
                     'T1_DcA','T2_DfB','T3_Coupling','DenomSum', ...
                     'T1_pct','T2_pct','T3_pct'}));
fprintf('  (T1_pct / T2_pct / T3_pct = each term as %% of the denominator sum)\n');
fprintf('--------------------------------------------------------------------------\n');

%% 6. PLOT
figure('Color','w','Position',[100 100 560 540]); hold on; box on;
allv=[exp_lives pred_int pred_lin]; allv=allv(allv>0);
lim=[10^floor(log10(min(allv))) 10^ceil(log10(max(allv)))];
plot(lim,lim,'k-','LineWidth',1.4); plot(lim,2*lim,'k--'); plot(lim,0.5*lim,'k--');
plot(lim,3*lim,'b-.'); plot(lim,lim/3,'b-.');
hI=scatter(exp_lives,pred_int,90,'filled','MarkerFaceColor',[0.10 0.45 0.80],'MarkerEdgeColor','k');
hL=scatter(exp_lives,pred_lin,90,'d','MarkerEdgeColor','r','LineWidth',1.4);
for k=1:N, text(exp_lives(k)*1.07,pred_int(k),sprintf('%.3f',cases(k).epsilon_a),'FontSize',9); end
set(gca,'XScale','log','YScale','log','FontSize',12); xlim(lim); ylim(lim); axis square; grid on;
xlabel('Experimental life (cycles)'); ylabel('Predicted life (cycles)');
title('30/30 constants on 0/60 data');
legend([hI hL],{'Interaction (imported)','Linear rule'},'Location','northwest');
saveas(gcf,'robustness_0_60.png');

%% ===================== HELPERS ==========================================
function [Dc,dSf]=damage(t,sv,B1,n1,Q,R,T,cw,fw)
    phi=B1*exp(-Q/(R*T)); Dc=0; dSf=0;
    for i=2:numel(t)
        dt=t(i)-t(i-1); s=sv(i); tc=t(i);
        if inwin(tc,cw)&&s>0, Dc=Dc+(1/phi)*(s)^(1-n1)*dt;
        elseif inwin(tc,fw),  dSf=dSf+s*dt; end
    end
end
function tf=inwin(tc,W)
    tf=false;
    for r=1:size(W,1), if tc>W(r,1)&&tc<=W(r,2), tf=true; return; end, end
end