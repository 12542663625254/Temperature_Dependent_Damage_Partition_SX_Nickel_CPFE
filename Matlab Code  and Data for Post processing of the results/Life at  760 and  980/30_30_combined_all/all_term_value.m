clc; clear; close all;
% =========================================================================
%  STAGE 2 - FIT a, d, C, Qd ON 30/30 DWELL DATA  (b fixed from Stage 1)
%
%     N_cf = 1/( Dc^a + Df^b + M*(Dc*Df)^d ),   M = C*exp(-Qd/RT)
%
%  C and Qd are separable ONLY if the 30/30 data spans >=2 temperatures.
%  Put the 1033 K cases AND the 1253 K (980 C) cases below.
% =========================================================================

%% ---- CAPTURE EVERYTHING PRINTED TO THE COMMAND WINDOW ----
diary_txt = 'stage2_command_window.txt';
diary_pdf = 'stage2_command_window.pdf';
if isfile(diary_txt), delete(diary_txt); end
diary(diary_txt); diary on;

%% ---- FROM STAGE 1 ----
b = 0.8807;          % fatigue exponent fitted on pure-fatigue data (FIXED here)

R  = 8.314;          % from table
n1 = 0.6;  Q = 6.97e-19;
Sg = 0.38; Sc = 0.9*Sg; Dfc = 0.9; D0 = 0;

%% 1. 30/30 CASES  (each tagged with its temperature in KELVIN)
% mk(file, ea, exp_life, T_K, creepWindows[2x2], fatigueWindows[2x2])
c = struct('file',{},'ea',{},'exp',{},'T',{},'cw',{},'fw',{});

% ---- 760 C = 1033 K (your 5 cases) ----
c(1)=mk('shi_30_30_2_4.csv',0.012,5,   1033,[12 42;66 96],   [42 66;96 120]);
c(2)=mk('pin_30_30_2_0.csv',0.010,210, 1033,[10 40;60 90],   [40 60;90 110]);
c(3)=mk('Fatigue.csv',      0.009,370, 1033,[681 711;729 759],[711 729;759 777]);
c(4)=mk('shi_30_30_1_6.csv',0.008,562, 1033,[560 590;606 636],[590 606;636 652]);
c(5)=mk('30_30_0_007.csv',  0.007,1455,1033,[7 37;51 81],     [37 51;81 95]);

% ---- 980 C = 1253 K : FILL IN the four 30/30 cases (lives 10,65,220,465) ----
% c(6)=mk('<file>',<ea>,10, 1253,[<cw1a cw1b; cw2a cw2b>],[<fw...>]);
% c(7)=mk('<file>',<ea>,65, 1253,[...],[...]);
% c(8)=mk('<file>',<ea>,220,1253,[...],[...]);
% c(9)=mk('<file>',<ea>,465,1253,[...],[...]);

N=numel(c);

%% 2. COMPUTE Dc, Df FOR EVERY CASE
Dc=zeros(1,N); Df=zeros(1,N); Tvec=zeros(1,N); exp_lives=zeros(1,N); ea_vec=zeros(1,N);
for k=1:N
    if ~isfile(c(k).file), error('File %s not found!',c(k).file); end
    raw=readmatrix(c(k).file);
    B1 = 2.04e3 - 2.11e4*c(k).ea - 1.23e7*c(k).ea^2;   % table formula
    [Dc(k),dSf]=damage(raw(:,1),raw(:,25),B1,n1,Q,R,c(k).T,c(k).cw,c(k).fw);
    if dSf>=Sg, Df(k)=1.0;
    else, Df(k)=max(0,D0+((Dfc-D0)/log(1-Sc/Sg))*log(max(1e-10,1-dSf/Sg))); end
    Tvec(k)=c(k).T; exp_lives(k)=c(k).exp; ea_vec(k)=c(k).ea;
end
disp(table((1:N)',Tvec',ea_vec',Dc',Df',exp_lives', ...
    'VariableNames',{'Case','T_K','Strain','Dc','Df','Exp_Life'}));

multiT = numel(unique(Tvec))>=2;
if ~multiT
    warning(['Only ONE temperature present -> C and Qd are NOT separable. ' ...
             'The fit will report M only; add the 1253 K cases to get C and Qd.']);
end

%% 3. FIT [a, d, C, Qd]  (b fixed) WITH MULTI-START
opt=optimset('Display','off','MaxFunEvals',2e4,'MaxIter',2e4,'TolX',1e-9,'TolFun',1e-12);
best_err=inf; best_p=[]; rng(1);
for s=1:120
    x0=[0.5+1.5*rand, 0.5+1.5*rand, 10^(4+3*rand), 5e5*rand];   % [a d C Qd]
    [p,e]=fminsearch(@(p) gcost(p,b,Dc,Df,Tvec,R,exp_lives),x0,opt);
    if e<best_err, best_err=e; best_p=p; end
end
a=best_p(1); d=best_p(2); C=best_p(3); Qd=best_p(4);
M = C*exp(-Qd./(R*Tvec));
pred = 1 ./ (Dc.^a + Df.^b + M.*(Dc.*Df).^d);
err  = abs(pred-exp_lives)./exp_lives*100;

%% 3b. BREAK DOWN THE THREE DENOMINATOR TERMS PER CASE
term1 = Dc.^a;                 % D_c^p  (creep term)
term2 = Df.^b;                 % D_f^q  (fatigue term)
term3 = M.*(Dc.*Df).^d;        % M*(Dc*Df)^d  (interaction term)
denom = term1 + term2 + term3; % = 1/pred

frac1 = term1./denom*100;      % % contribution of each term
frac2 = term2./denom*100;
frac3 = term3./denom*100;

%% 4. REPORT
fprintf('\n=============== STAGE 2: 30/30 CONSTANTS ===============\n');
fprintf(' a  = %.4f\n b  = %.4f  (fixed from Stage 1)\n d  = %.4f\n', a,b,d);
if multiT
    fprintf(' C  = %.6g\n Qd = %.1f J/mol (= %.3f kJ/mol)\n', C, Qd, Qd/1000);
else
    fprintf(' M  = %.6g  (lumped; C,Qd not separable at one temperature)\n', C*exp(-Qd/(R*Tvec(1))));
end
fprintf(' global log10-RMSE = %.4f\n', best_err);
fprintf('========================================================\n');
disp(table((1:N)',Tvec',ea_vec',exp_lives',pred',err', ...
    'VariableNames',{'Case','T_K','Strain','Exp_Life','Pred_Life','Error_Pct'}));

%% 4b. DENOMINATOR TERM BREAKDOWN
fprintf('\n=========== DENOMINATOR TERM BREAKDOWN ===========\n');
fprintf('  N_cf = 1 / ( Dc^p + Df^q + M*(Dc*Df)^d )\n\n');
format shortG
disp(table((1:N)', Tvec', ea_vec', term1', term2', term3', denom', ...
    'VariableNames',{'Case','T_K','Strain','Dc_p','Df_q','Interaction','Denom_Total'}));

fprintf('\n--- Percentage contribution of each term ---\n');
disp(table((1:N)', Tvec', ea_vec', frac1', frac2', frac3', ...
    'VariableNames',{'Case','T_K','Strain','Creep_pct','Fatigue_pct','Interaction_pct'}));
format short

%% PLOT
figure('Color','w','Position',[100 100 560 540]); hold on; box on;
allv=[exp_lives pred]; allv=allv(allv>0);
lim=[10^floor(log10(min(allv))) 10^ceil(log10(max(allv)))];
plot(lim,lim,'k-','LineWidth',1.4); plot(lim,2*lim,'k--'); plot(lim,0.5*lim,'k--');
mk1033=Tvec==1033;
scatter(exp_lives(mk1033),pred(mk1033),90,'filled','MarkerFaceColor',[0.10 0.45 0.80],'MarkerEdgeColor','k');
scatter(exp_lives(~mk1033),pred(~mk1033),90,'filled','MarkerFaceColor',[0.85 0.33 0.10],'MarkerEdgeColor','k');
for k=1:N, text(exp_lives(k)*1.07,pred(k),sprintf('%.3f',ea_vec(k)),'FontSize',9); end
set(gca,'XScale','log','YScale','log','FontSize',12); xlim(lim); ylim(lim); axis square; grid on;
xlabel('Experimental life (cycles)'); ylabel('Predicted life (cycles)');
title('Stage 2: 30/30 fit'); legend({'Perfect','x2 / /2','','1033 K','1253 K'},'Location','northwest');
saveas(gcf,'stage2_30_30.png');

%% ---- STOP CAPTURING AND CONVERT THE LOG TO PDF ----
diary off;
save_diary_as_pdf(diary_txt, diary_pdf);

%% ===================== HELPERS ==========================================
function s=mk(file,ea,exp,T,cw,fw)
    s.file=file; s.ea=ea; s.exp=exp; s.T=T; s.cw=cw; s.fw=fw;
end

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

function err=gcost(p,b,Dc,Df,Tvec,R,exp_lives)
    a=p(1); d=p(2); C=p(3); Qd=p(4); pen=0;
    if a<0.3||a>3, pen=pen+1e7; end
    if d<0.3||d>3, pen=pen+1e7; end
    if C<=0,       pen=pen+1e7; end
    if Qd<0||Qd>6e5, pen=pen+1e7; end
    if pen>0, err=1e10+pen; return; end
    M=C*exp(-Qd./(R*Tvec));
    pl=1 ./ (Dc.^a + Df.^b + M.*(Dc.*Df).^d);
    err=sqrt(mean((log10(exp_lives)-log10(max(1e-12,pl))).^2));
end

% --- Converts a plain-text diary log into a paginated PDF -----------------
% Renders the text into monospaced figure pages (Courier) and writes them
% to one multi-page PDF using print(...,'-append'). Works on any platform
% with no extra toolboxes.
function save_diary_as_pdf(diary_file, pdf_file, lines_per_page)
    if nargin < 3, lines_per_page = 65; end
    if ~isfile(diary_file)
        warning('Diary file %s not found, skipping PDF export.', diary_file);
        return;
    end
    txt = fileread(diary_file);
    txt = strrep(txt, char(13), '');           % strip CR if present
    rawlines = strsplit(txt, char(10));
    nLines = numel(rawlines);
    nPages = max(1, ceil(nLines/lines_per_page));
    if isfile(pdf_file), delete(pdf_file); end
    for p = 1:nPages
        idx0 = (p-1)*lines_per_page + 1;
        idx1 = min(p*lines_per_page, nLines);
        pageLines = rawlines(idx0:idx1);
        pageStr = strjoin(pageLines, newline);
        fig = figure('Color','w','Visible','off','Units','inches', ...
                      'Position',[0 0 8.5 11],'PaperUnits','inches', ...
                      'PaperSize',[8.5 11],'PaperPosition',[0 0 8.5 11]);
        axes('Parent',fig,'Position',[0 0 1 1],'Visible','off'); hold on;
        text(0.02,0.99,pageStr,'Units','normalized', ...
             'VerticalAlignment','top','HorizontalAlignment','left', ...
             'FontName','Courier New','FontSize',7.5,'Interpreter','none');
        if p == 1
            print(fig, pdf_file, '-dpdf', '-r150');
        else
            print(fig, pdf_file, '-dpdf', '-r150', '-append');
        end
        close(fig);
    end
    fprintf('Saved command-window output to %s (%d page(s)).\n', pdf_file, nPages);
end