function [ IntVar ] = update_aging_metrics_3Para(Config, IntVar)

%[IntVar.Swell_alpha, IntVar.Swell_beta, IntVar.Cap_fade_alpha, IntVar.Cap_fade_beta, IntVar.WRc_alpha, IntVar.WRc_beta, IntVar.WRa_alpha, IntVar.WRa_beta] = Get_Aging_Param_T_now_v4(Config.opt_para_all,IntVar.T_all);
Vcath_now = IntVar.Vcath_all;
Vref_now = IntVar.Vref_all;
t_Vref = IntVar.t_all;

%% Hard swell

% VV = (ones(size(Vcath_now))*Config.opt_para_all.swell.alpha_all').*exp(-(Vcath_now - Config.V_th)*Config.opt_para_all.swell.beta_all')+ ...
%      ((IntVar.I_all>0).*IntVar.I_all.*exp(Vref_now*Config.beta2_fixed))*Config.opt_para_all.swell.gamma_all';
% [R_swell] = eval_aging_rate(IntVar.T_all,Config.opt_para_all.swell.Temp_grid, VV);
% % min(VV(:,2)-VV(:,1))
% % min(VV(:,3)-VV(:,2)) 
% aging_rate_this_cycle_hardswell = R_swell.*(IntVar.t_clock + t_Vref).^(-0.5);
% d_hard_swell_this_cycle = trapz(t_Vref,aging_rate_this_cycle_hardswell);
% IntVar.d_hard_swell_now = IntVar.d_hard_swell_now + d_hard_swell_this_cycle;
%% Cap Fade

VV = (ones(size(Vcath_now))*Config.opt_para.cap_fade.alpha_all').*exp(-(Vcath_now - Config.V_th)*Config.opt_para.cap_fade.beta1_all')+ ...
     ((IntVar.I_all>0).*IntVar.I_all.*exp(Vref_now*Config.beta2_fixed))*Config.opt_para_all.cap_fade.gamma_all';
% min(VV(:,1)-VV(:,2))
% min(VV(:,2)-VV(:,3)) 

[R_CF] = eval_aging_rate(IntVar.T_all,Config.opt_para_all.cap_fade.Temp_grid, VV);

CF_rate_this_cyle = trapz(t_Vref,R_CF);

IntVar.cap_fade_now = IntVar.cap_fade_now + CF_rate_this_cyle;
if abs(IntVar.cap_fade_now)>50
   disp('Cap. Fade is larger than 50%. Is it correct?') 
end
IntVar.Cap_now = Config.Cap_new*(100+IntVar.cap_fade_now)/100;

%% WRc growth
VV = (ones(size(Vcath_now))*Config.opt_para_all.WRc_growth.alpha_all').*exp(-(Vcath_now - Config.V_th)*Config.opt_para_all.WRc_growth.beta_all')+ ...
     ((IntVar.I_all>0).*IntVar.I_all.*exp(Vref_now*Config.beta2_fixed))*Config.opt_para_all.WRc_growth.gamma_all';
[R_WRc] = eval_aging_rate(IntVar.T_all, Config.opt_para_all.WRc_growth.Temp_grid, VV);

% min(VV(:,2)-VV(:,1))
% min(VV(:,3)-VV(:,2)) 

WRc_growth_rate_this_cyle = trapz(t_Vref,R_WRc);
IntVar.dWRc_now_in_percent = IntVar.dWRc_now_in_percent+WRc_growth_rate_this_cyle;

%% WRa growth
VV = (ones(size(Vcath_now))*Config.opt_para_all.WRa_growth.alpha_all').*exp(-(Vcath_now - Config.V_th)*Config.opt_para_all.WRa_growth.beta_all')+ ...
     ((IntVar.I_all>0).*IntVar.I_all.*exp(Vref_now*Config.beta2_fixed))*Config.opt_para_all.WRa_growth.gamma_all';
[R_WRa] = eval_aging_rate(IntVar.T_all, Config.opt_para_all.WRa_growth.Temp_grid, VV);

% min(VV(:,2)-VV(:,1))
% min(VV(:,3)-VV(:,2)) 
    
WRa_growth_rate_this_cyle = trapz(t_Vref,R_WRa);
IntVar.dWRa_now_in_percent = IntVar.dWRa_now_in_percent+WRa_growth_rate_this_cyle;

%% gaseous swell model
% IntVar = sim_gas_swell_model_RDC(IntVar,Config);

% figure(300001)
% plot(t_Vref/3600, R_WRa)
% hold on
% indx = find(IntVar.I_all>0);
% figure(300001)
% plot(Vref_now(indx), IntVar.I_all(indx))
% hold on
%%
IntVar.t_clock = IntVar.t_clock + t_Vref(end) - t_Vref(1);
end