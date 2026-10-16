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
data.tWw_f = [ ... % time since metam (d), wet weight (g)
292	0.53
292	0.66
292	0.54
292	0.65
292	0.47
292	0.57
292	0.67
292	0.61
292	0.52
292	0.47
292	0.69
292	0.58
292	0.44
292	0.7
292	0.54
292	0.49
292	0.86
292	0.53
658	1.196
658	1.844
658	1.201
658	1.462
658	1.73
658	1.214
658	1.458
658	1.13
658	1.118
658	1.202
658	1.783
658	1.545
658	0.822
658	1.813
658	1.275
658	1.157
658	2.042
658	1.119
1007	2.939
1008	3.561
1007	3.241
1006	3.348
1006	3.321
1006	3.304
1006	2.878
1007	2.596
1007	2.779
1006	3.461
1006	3.569
1006	2.941
1007	4.036
1008	2.899
1008	1.8
1007	3.367
1007	3.574
1315	4.28
1316	3.76
1316	3.75
1316	3.95
1315	3.51
1316	3.88
1315	4.32
1315	3.21
1316	3.05
1316	3.11
1317	4.79
1316	2.28
1316	3.62
1317	2.72
1315	2.37
1315	2.72
1315	4.48
1315	2.69
1663	4.06
1663	4.61
1663	4.47
1663	4.13
1663	4.45
1663	4.66
1664	4.28
1664	3.76
1664	3.59
1664	4.26
1664	4.99
1664	4.35
1664	3.93
1664	3.57
1665	4.49
1665	4.32
1665	4.88];
for i=2:size(data.tWw_f,1); if data.tWw_f(i,1)<=data.tWw_f(i-1,1); data.tWw_f(i,1)=data.tWw_f(i-1,1)+1e-3; end; end
units.tWw_f = {'d', 'g'}; label.tWw_f = {'time since metam', 'wet weight', 'females'};
temp.tWw_f = C2K(12); units.temp.tWw_f = 'K'; label.temp.tWw_f = 'temperature';
comment.tWw_f = 'Data for females'; bibkey.tWw_f = 'BaskGvoz2021';
%
data.tWw_m = [ ... % time since metam (d), wet weight (g)
292	0.73
292	0.63
292	0.48
292	0.5
292	0.64
292	0.52
292	0.58
292	0.62
292	0.53
292	0.51
292	0.52
292	0.63
292	0.65
292	0.4
292	0.61
292	0.65
292	0.49
292	0.83
292	0.46
658	2.032
658	1.577
658	1.452
658	1.59
658	1.47
658	1.571
658	1.529
658	1.705
658	1.055
658	1.21
658	1.229
658	1.546
658	1.679
658	0.885
658	1.516
658	1.472
658	1.88
658	1.282
1007	3.004
1007	3.978
1007	2.754
1008	2.601
1006	2.811
1006	2.145
1006	2.389
1007	2.408
1006	2.215
1007	2.358
1006	2.135
1007	2.522
1006	2.449
1006	1.956
1007	2.446
1008	2.202
1008	2.24
1007	2.588
1007	1.8
1316	2.64
1316	4.94
1316	2.85
1315	3.34
1316	2.6
1316	2.7
1315	2.47
1316	2.5
1315	2.45
1316	2.17
1315	2.55
1317	2.38
1317	3.92
1315	3.23
1317	3.15
1315	3.95
1315	2.17
1663	2.74
1663	2.75
1663	2.71
1663	3.41
1663	2.74
1663	2.91
1663	2.71
1663	2.63
1663	2.65
1664	2.11
1664	2.34
1664	2.53
1664	2.5
1664	2.81
1664	2.61
1665	2.86
1665	2.76];
for i=2:size(data.tWw_m,1); if data.tWw_m(i,1)<=data.tWw_m(i-1,1); data.tWw_m(i,1)=data.tWw_m(i-1,1)+1e-3; end; end
units.tWw_m = {'d', 'g'}; label.tWw_m = {'time since metam', 'wet weight' 'males'};
temp.tWw_m = C2K(12); units.temp.tWw_m = 'K'; label.temp.tWw_m = 'temperature';
comment.tWw_m = 'Data for males'; bibkey.tWw_m = 'BaskGvoz2021';

% time-respiration
data.tJO_f = [ ... % time since birth (d), O2 consumption (mL/h)
292             0.11495514
292             0.10307827
292             0.08463952
292             0.11467741
292             0.04027615
292             0.09804198
292             0.12709846
292             0.05067992
292             0.04428620
292             0.10506425
292             0.12528652
292             0.07584610
292             0.04778457
292             0.10935149
292             0.11750344
292             0.04264609
292             0.15991204
292             0.05591683
658             0.16750800
658             0.09941124
658             0.06118746
658             0.21047326
658             0.06302341
658             0.12780474
658             0.21783045
658             0.06488677
658             0.12383814
658             0.07430728
658             0.11122167
658             0.17464138
658             0.05771578
658             0.18815208
658             0.14520316
658             0.05492558
658             0.20104789
658             0.08880951
1007             0.40386940
1008             0.16548833
1007             0.40406752
1006             0.24494644
1006             0.21199208
1006             0.38177561
1006             0.27093334
1007             0.30574163
1007             0.25565725
1006             0.20717700
1006             0.20320797
1006             0.21186675
1007             0.40249458
1008             0.34332529
1008             0.17425308
1007             0.27945370
1007             0.22189487
1315             0.20574450
1316             0.41278664
1316             0.36723752
1316             0.40108717
1315             0.20097966
1316             0.20258585
1315             0.48661086
1315             0.17916698
1316             0.39608019
1316             0.46345937
1317             0.45596884
1316             0.19972721
1316             0.21659993
1317             0.15628415
1315             0.29088750
1315             0.15333433
1315             0.35969210
1315             0.38208050
1663             0.43325024
1663             0.34188272
1663             0.41861139
1663             0.29729225
1663             0.25355986
1663             0.26442675
1664             0.23790704
1664             0.25643493
1664             0.30882147
1664             0.18140554
1664             0.39878737
1664             0.22916219
1664             0.21439196
1664             0.48932641
1665             0.23686263
1665             0.39573906
1665             0.52202144];
for i=2:size(data.tJO_f,1); if data.tJO_f(i,1)<=data.tJO_f(i-1,1); data.tJO_f(i,1)=data.tJO_f(i-1,1)+1e-3; end; end
units.tJO_f  = {'d', 'mL/h'};  label.tJO_f = {'time since birth', 'O2 consumption','females'};  
temp.tJO_f   = C2K(12);  units.temp.tJO_f = 'K'; label.temp.tJO_f = 'temperature';
bibkey.tJO_f = 'BaskGvoz2021';
comment.tJO_f = 'data for females';
%
data.tJO_m = [ ... % time since birth (d), O2 consumption (mL/h)
292             0.13367042
292             0.05646998
292             0.05298743
292             0.08601500
292             0.11893019
292             0.07219962
292             0.05387069
292             0.11232635
292             0.08165338
292             0.11679215
292             0.10442352
292             0.06630379
292             0.14164141
292             0.04339970
292             0.12552492
292             0.09805789
292             0.13414897
292             0.16117562
292             0.07436749
658             0.09938976
658             0.08090096
658             0.16961307
658             0.22350147
658             0.13621321
658             0.17447202
658             0.08881616
658             0.13277801
658             0.13170438
658             0.23073877
658             0.18527416
658             0.09292764
658             0.26098328
658             0.13508583
658             0.18580809
658             0.08671667
658             0.13085192
658             0.19714807
1007             0.19458054
1007             0.36450933
1007             0.42655797
1008             0.07770985
1006             0.14786943
1006             0.26589356
1006             0.18245941
1007             0.13087229
1006             0.24342475
1007             0.37126120
1006             0.31359709
1007             0.12940252
1006             0.40361661
1006             0.35385191
1007             0.27140033
1008             0.12600380
1008             0.20612697
1007             0.37515752
1007             0.34032693
1316             0.21789756
1316             0.28790995
1316             0.34984560
1315             0.41034468
1316             0.30846511
1316             0.16257203
1315             0.27393706
1316             0.17325685
1315             0.32600176
1316             0.20556823
1315             0.37671202
1317             0.20786018
1317             0.16624099
1315             0.46511786
1317             0.16259450
1315             0.40268991
1315             0.13484217
1663             0.29902492
1663             0.42798275
1663             0.16757513
1663             0.43814959
1663             0.30753906
1663             0.16544639
1663             0.29081439
1663             0.33036958
1663             0.20735054
1664             0.14272896
1664             0.39004210
1664             0.26675571
1664             0.16303780
1664             0.35135893
1664             0.29933925
1665             0.15482916
1665             0.16219263];
for i=2:size(data.tJO_m,1); if data.tJO_m(i,1)<=data.tJO_m(i-1,1); data.tJO_m(i,1)=data.tJO_m(i-1,1)+1e-3; end; end
units.tJO_m  = {'d', 'mL/h'};  label.tJO_m = {'time since birth', 'O2 consumption','males'};  
temp.tJO_m   = C2K(12);  units.temp.tJO_m = 'K'; label.temp.tJO_m = 'temperature';
bibkey.tJO_m = 'BaskGvoz2021';
comment.tJO_m = 'data for males';

%% set weights for all real data
weights = setweights(data, []);
weights.tL_f = 5 * weights.tL_f;
weights.tL_m = 5 * weights.tL_m;
weights.tWw_f = 5 * weights.tWw_f;
weights.tWw_m = 5 * weights.tWw_m;

%% set pseudodata and respective weights
[data, units, label, weights] = addpseudodata(data, units, label, weights);

%% pack auxData and txtData for output
auxData.temp = temp;
txtData.units = units;
txtData.label = label;
txtData.bibkey = bibkey;
txtData.comment = comment;

%% Group plots
set1 = {'tL_f','tL_m'}; subtitle1 = {'Data for females, males'};
set2 = {'tWw_f','tWw_m'}; subtitle2 = {'Data for females, males'};
set3 = {'tJO_f','tJO_m'}; subtitle3 = {'Data for females, males'};
metaData.grp.sets = {set1, set2, set3};
metaData.grp.subtitle = {subtitle1, subtitle2, subtitle3};

%% Discussion points
D1 = 'Males are assumed to differ from females by {p_Am} and E_Hp only';
D2 = 'Temperatures are guessed';
D3 = 'mod_1: males have equal state variables at b, compared to females';
D4 = 'mod_2: tWw and tJO data added';
metaData.discussion = struct('D1',D1, 'D2',D2, 'D3',D3, 'D4',D4);  

%% Facts
% F1 = '';
% metaData.bibkey.F1 = ''; 
% metaData.facts = struct('F1',F1);

%% Links
metaData.links.id_CoL = '6N65G'; % Cat of Life
metaData.links.id_ITIS = '775917'; % ITIS
metaData.links.id_EoL = '330758'; % Ency of Life
metaData.links.id_Wiki = 'Ichthyosaura_alpestris'; % Wikipedia
metaData.links.id_ADW = 'Ichthyosaura_alpestris'; % ADW
metaData.links.id_Taxo = '1397948'; % Taxonomicon
metaData.links.id_WoRMS = '1663441'; % WoRMS
metaData.links.id_amphweb = 'Ichthyosaura+alpestris'; % AmphibiaWeb
metaData.links.id_AnAge = 'Ichthyosaura_alpestris'; % AnAge

%% References
bibkey = 'Wiki'; type = 'Misc'; bib = ...
'howpublished = {\url{https://en.wikipedia.org/wiki/Ichthyosaura_alpestris}}';
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Kooy2010'; type = 'Book'; bib = [ ...  % used in setting of chemical parameters and pseudodata
'author = {Kooijman, S.A.L.M.}, ' ...
'year = {2010}, ' ...
'title  = {Dynamic Energy Budget theory for metabolic organisation}, ' ...
'publisher = {Cambridge Univ. Press, Cambridge}, ' ...
'pages = {Table 4.2 (page 150), 8.1 (page 300)}, ' ...
'howpublished = {\url{../../../bib/Kooy2010.html}}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'MiauGuye2000'; type = 'article'; bib = [ ...
'author = {Claude Miaud and Robert Guyetant and Helmut Faber}, ' ... 
'year = {2000}, ' ...
'title = {Age, size, and growth of the {A}lpine newt, \emph{Triturus alpestris} ({U}rodela: {S}alamandridae), at high altitude and a review of life-history trait variation throughout its range}, ' ...
'journal = {Herpetologica}, ' ...
'volume = {56(2)}, ' ...
'pages = {135-144}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'BaskGvoz2021'; type = 'article'; bib = [ ...
'author = {Senka Ba\v{s}kiera and Lum\''{i}r Gvo\v{z}d\''{i}k}, ' ... 
'year = {2021}, ' ...
'title = {Repeatability and heritability of resting metabolic rate in a long-lived amphibian}, ' ...
'journal = {Comparative Biochemistry and Physiology Part A: Molecular \& Integrative Physiology}, ' ...
'volume = {253}, ' ...
'pages = {110858}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Gvoz2026'; type = 'article'; bib = [ ...
'author = {Lumír Gvoždík}, ' ... 
'title = {Unpublished data}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
