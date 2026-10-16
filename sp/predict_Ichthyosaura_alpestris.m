function [prdData, info] = predict_Ichthyosaura_alpestris(par, data, auxData)
  
  % unpack par, data, auxData
  cPar = parscomp_st(par); vars_pull(par); 
  vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);


  if E_Hpm < E_Hb
    prdData = []; info = 0; return
  end

  % compute temperature correction factors
  TC = tempcorr(temp.am, T_ref, T_A); kT_M = TC * k_M;
  TC_tL = tempcorr(temp.tL_f, T_ref, T_A); 
  TC_tW = tempcorr(temp.tWw_f, T_ref, T_A);

  % zero-variate data

  % life cycle
  pars_tj = [g; k; l_T; v_Hb; v_Hj];               % compose parameter vector
  [tau_j, tau_b, l_j, l_b] = get_tp(pars_tj, f);   % -, scaled times & lengths at f
  pars_tp = [g; k; l_T; v_Hb; v_Hp];               % compose parameter vector
  [tau_p, tau_b, l_p, l_b, info] = get_tp(pars_tp, f); % -, scaled times & lengths at f
  
  % birth
  L_b = L_m * l_b;                  % cm, structural length at birth at f
  Ww_b = L_b^3 * (1 + f * w);       % g, wet weight at birth at f
  
  % metam
  L_j = L_m * l_j;                  % cm, struc length at metam
  Lw_j = L_j/ del_M;                % cm, SVL at metam

  % puberty 
  L_p = L_m * l_p;                  % cm, struc length at puberty
  Lw_p = L_p/ del_M;                % cm, SVL at puberty
  tT_p = (tau_p - tau_b)/ kT_M;     % d, time since birth at puberty

  % ultimate
  l_i = f - l_T;                    % -, scaled ultimate length
  L_i = L_m * l_i;                  % cm, ultimate structural length at f
  Lw_i = L_i/ del_M;                % cm, ultimate SVL at f
  Ww_i = L_i^3 * (1 + f * w);       % g, ultimate wet weight

  % reproduction
  pars_R = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hp]; % compose parameter vector at T
  RT_i = TC *  reprod_rate(L_i, f, pars_R);% #/d, reprod rate 

  % life span
  pars_tm = [g; l_T; h_a/ k_M^2; s_G];  % compose parameter vector at T_ref
  t_m = get_tm_s(pars_tm, f, l_b);      % -, scaled mean life span at T_ref
  aT_m = t_m/ kT_M;                     % d, mean life span at T
  
  % males
  p_Am_m = z_m * p_M/ kap;             % J/d.cm^2, {p_Am} spec assimilation flux
  E_m_m = p_Am_m/ v;                   % J/cm^3, reserve capacity [E_m]
  g_m = E_G/ (kap* E_m_m);             % -, energy investment ratio
  m_Em_m = y_E_V * E_m_m/ E_G;         % mol/mol, reserve capacity 
  ome_m = m_Em_m * w_E/ w_V;           % -, contribution of reserve to weight
  L_mm = v/ k_M/ g_m;                  % cm, max struct length
  L_im = f * L_mm; Lw_im = L_im/ del_M;% cm, ultimate struct, SVL
  Ww_im = L_im^3 * (1 + f * ome_m);    % g, ultimate wet weight
  pars_tpm = [g_m k l_T v_Hb v_Hpm]; pars_tjm = [g_m k l_T v_Hb v_Hj];
  [t_pm, t_bm, l_pm, l_bm] = get_tp(pars_tpm, f);
  L_pm = L_mm * l_pm; Lw_pm = L_pm/ del_M; % cm, SVL at puberty
  tT_pm = (t_pm - t_bm)/ kT_M;   % d, time since birth at puberty for males
  Lw_im = L_im/ del_M;      % cm, ultimate SVL

  
  % pack to output
  prdData.tp = tT_p;
  prdData.tpm = tT_pm;
  prdData.am = aT_m;
  prdData.Lj = Lw_j;
  prdData.Lp = Lw_p;
  prdData.Lpm = Lw_pm;
  prdData.Li = Lw_i;
  prdData.Lim = Lw_im;
  prdData.Wwb = Ww_b;
  prdData.Wwi = Ww_i;
  prdData.Wwim = Ww_im;
  prdData.Ri = RT_i;
  
  % univar data
  
  % time since metam-length
  % female
  tvel = get_tp(pars_tp, f_tL, get_lp(pars_tj, f_tL), tL_f(:,1)*k_M*TC_tL);
  ELw_f = L_m * tvel(:,4)/ del_M; % cm, SVL length 
  % male
  tvel = get_tp(pars_tpm, f_tL, get_lp(pars_tjm, f_tL), tL_m(:,1)*k_M*TC_tL);
  ELw_m = L_mm * tvel(:,4)/ del_M; % cm, SVL length 

  % time since metam-weight
  % female
  tvel = get_tp(pars_tp, f_tW, get_lp(pars_tj, f_tW), tWw_f(:,1)*k_M*TC_tW);
  EWw_f = (L_m * tvel(:,4)).^3 .* (1 + ome * tvel(:,3)); % g, weight 
  % male
  tvel = get_tp(pars_tpm, f_tW, get_lp(pars_tjm, f_tW), tWw_m(:,1)*k_M*TC_tW);
  EWw_m = (L_mm * tvel(:,4)).^3 .* (1 + ome_m * tvel(:,3)); % g, weight 
  
  
  % time-respiration
  p_ref = p_Am * L_m^2;             % J/d, max assimilation power at max size
  pars_power = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hp];
  tvel = get_tp(pars_tp, f_tW, [], tJO_f(:,1)*k_M*TC_tW);
  L = L_m * tvel(:,4);  % cm, struct length
  pACSJGRD = p_ref * scaled_power(L, f_tW, pars_power, l_b, l_p);
  pADG = pACSJGRD(:, [1 7 5])'; pADG(1,:) = 0;  % J/d, select assim, dissip, growth, set assim to zero
  J_O = eta_O * pADG;               % J_X, J_V, J_E, J_P in rows, A, D, G in cols
  J_M = - (n_M\n_O) * J_O;          % J_C, J_H, J_O, J_N in rows, A, D, G in cols
  X_gas = 1/ 24.4;     % M, mol of gas per litre at 20 C and 1 bar 
  EJO_f = - 1e3 * J_M(3,:)'/ 24/ X_gas;  % ml/h, dioxygen use
  %
  p_refm = p_Am_m * L_mm^2;         % J/d, max assimilation power at max size
  pars_power_m = [kap; kap_R; g_m; k_J; k_M; L_T; v; U_Hb; U_Hp];
  tvel = get_tp(pars_tpm, f_tW, [], tJO_m(:,1)*k_M*TC_tW);
  L = L_mm * tvel(:,4);  % cm, struct length
  pACSJGRD = p_refm * scaled_power(L, f_tW, pars_power_m, l_bm, l_pm);
  pADG = pACSJGRD(:, [1 7 5])'; pADG(1,:) = 0;  % J/d, select assim, dissip, growth, set assim to zero
  J_O = eta_O * pADG;               % mol/d, J_X, J_V, J_E, J_P in rows, A, D, G in cols
  J_M = - (n_M\n_O) * J_O;          % mol/d, J_C, J_H, J_O, J_N in rows, A, D, G in cols
  X_gas = 1/ 24.4;     % M, mol of gas per litre at 20 C and 1 bar 
  EJO_m = - 1e3 * J_M(3,:)'/ 24/ X_gas;  % ml/h, dioxygen use
  
  % pack to output
  prdData.tL_f = ELw_f;
  prdData.tL_m = ELw_m;
  prdData.tWw_f = EWw_f;
  prdData.tWw_m = EWw_m;
  prdData.tJO_f = EJO_f;
  prdData.tJO_m = EJO_m;
