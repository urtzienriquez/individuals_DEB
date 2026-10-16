function [data, auxData, metaData, txtData, weights] = mydata_Ichthyosaura_alpestris

%% set metaData
metaData.phylum     = 'Chordata'; 
metaData.class      = 'Amphibia'; 
metaData.order      = 'Caudata'; 
metaData.family     = 'Salamandridae';
metaData.species    = 'Ichthyosaura_alpestris'; 
metaData.species_en = 'Alpine newt'; 

metaData.ecoCode.climate = {'Cfb', 'Dfb'};
metaData.ecoCode.ecozone = {'THp'};
metaData.ecoCode.habitat = {'0jFp', 'jiTh', 'jiFm'};
metaData.ecoCode.embryo  = {'Fh'};
metaData.ecoCode.migrate = {'TT'};
metaData.ecoCode.food    = {'biCi'};
metaData.ecoCode.gender  = {'Dg'};
metaData.ecoCode.reprod  = {'O'};

metaData.T_typical  = C2K(12); % K, body temp 
metaData.data_0     = {'am'; 'Lj'; 'Lp'; 'Lpm'; 'Li'; 'Lim'; 'Wwb'; 'Wwi'; 'Wwim'; 'Ri'; 'tp'; 'tpm'}; 
metaData.data_1     = {'t-L'; 't-Ww'; 't-JO'}; 

metaData.COMPLETE = 3; % using criteria of LikaKear2011

metaData.author   = {'Bas Kooijman'};    
metaData.date_subm = [2019 06 12];              
metaData.email    = {'bas.kooijman@vu.nl'};            
metaData.address  = {'VU University, Amsterdam'};   

metaData.author_mod_1   = {'Bas Kooijman'};    
metaData.date_mod_1     = [2023 04 08];              
metaData.email_mod_1    = {'bas.kooijman@vu.nl'};            
metaData.address_mod_1  = {'VU University, Amsterdam'};   

metaData.author_mod_2   = {'Urtzi Enriquez-Urzelai'};    
metaData.date_mod_2     = [2026 07 20];              
metaData.email_mod_2    = {'uenriquez@proton.me'};            
metaData.address_mod_2  = {'Institute of Vertebrate Biology, Czech Academy of Sciences'};   

metaData.curator     = {'Bas Kooijman'};
metaData.email_cur   = {'bas.kooijman@vu.nl'}; 
metaData.date_acc    = [2026 07 21];

%% set data
% zero-variate data

data.am = 20 * 365;    units.am = 'd';    label.am = 'life span';                bibkey.am = 'MiauGuye2000';   
  temp.am = C2K(12);  units.temp.am = 'K'; label.temp.am = 'temperature'; 
data.tp  = 3 * 365;  units.tp  = 'd';  label.tp  = 'time since birth at puberty for females';  bibkey.tp  = 'MiauGuye2000';
  temp.tp  = C2K(12);  units.temp.tp  = 'K'; label.temp.tp  = 'temperature';
data.tpm = 2 * 365;  units.tpm = 'd';  label.tpm = 'time since birth at puberty for males';    bibkey.tpm = 'MiauGuye2000';
  temp.tpm = C2K(12);  units.temp.tpm = 'K'; label.temp.tpm = 'temperature';

data.Lj  = 1.83;    units.Lj  = 'cm';  label.Lj  = 'SVL at metam for females';   bibkey.Lj  = 'MiauGuye2000';
data.Lp  = 4.9;    units.Lp  = 'cm';  label.Lp  = 'SVL at puberty for females';   bibkey.Lp  = 'MiauGuye2000';
data.Lpm  = 4.4;    units.Lpm = 'cm';  label.Lpm = 'SVL at puberty for males';   bibkey.Lpm  = 'MiauGuye2000';
data.Li  = 5.84;    units.Li  = 'cm';  label.Li  = 'ultimate SVL for females';   bibkey.Li  = 'MiauGuye2000';
data.Lim = 5.05;    units.Lim = 'cm';  label.Lim = 'ultimate SVL for males';     bibkey.Lim  = 'MiauGuye2000';

data.Wwb = 3.5e-3; units.Wwb = 'g';   label.Wwb = 'wet weight at birth';        bibkey.Wwb = 'guess';
  comment.Wwb = 'based on Triturus cristatus';
data.Wwi = 4.65;    units.Wwi = 'g';   label.Wwi = 'ultimate wet weight for females';     bibkey.Wwi = 'Gvoz2026';
data.Wwim = 3.41;   units.Wwim = 'g';   label.Wwim = 'ultimate wet weight for males';     bibkey.Wwim = 'Gvoz2026';

data.Ri  = 600/365; units.Ri  = '#/d'; label.Ri  = 'maximum reprod rate';     bibkey.Ri  = 'Gvoz2026';   
  temp.Ri = C2K(12); units.temp.Ri = 'K'; label.temp.Ri = 'temperature';

% univariate data
% time-length
data.tL_f = [ ... % time since metam (yr), SVL (cm)
0.000	1.847
8.998	5.699
9.984	5.687
9.985	5.786
10.969	5.741
10.970	5.784
10.971	5.850
11.957	5.882
11.958	6.013
11.981	5.827
12.966	5.782
12.992	5.913
12.993	6.001
14.027	5.868
14.028	5.923
14.988	5.867
14.989	5.977
16.000	6.052
16.024	5.943
16.984	5.799
18.982	5.994];
data.tL_f(:,1) = data.tL_f(:,1) * 365; % convert yr to d 
units.tL_f  = {'d', 'cm'};  label.tL_f = {'time since metam', 'SVL', 'females'};  
temp.tL_f   = C2K(6);  units.temp.tL_f = 'K'; label.temp.tL_f = 'temperature';
bibkey.tL_f = 'MiauGuye2000';
comment.tL_f = 'Data for females';
%
data.tL_m = [ ... % time since metam (yr), SVL (cm)
0.000	1.891
11.011	4.898
11.012	4.953
11.038	5.139
11.972	4.941
11.998	5.039
11.999	5.171
12.984	5.071
13.008	5.038
13.009	5.137
13.994	5.059
13.995	5.146
15.004	4.992
15.992	5.199
16.015	5.089
16.953	5.187
19.960	5.282];
data.tL_m(:,1) = data.tL_m(:,1) * 365; % convert yr to d 
units.tL_m  = {'d', 'cm'};  label.tL_m = {'time since metam', 'SVL', 'males'};  
temp.tL_m   = C2K(6);  units.temp.tL_m = 'K'; label.temp.tL_m = 'temperature';
bibkey.tL_m = 'MiauGuye2000';
comment.tL_m = 'Data for males';

% time-weight
data.tWw = [ ... % time since metam (d), wet weight (g)
