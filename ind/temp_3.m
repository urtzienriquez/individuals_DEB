];
units.tJO  = {'d', 'mL/h'};  label.tJO = {'time since birth', 'O2 consumption'};
temp.tJO   = C2K(12);  units.temp.tJO = 'K'; label.temp.tJO = 'temperature';
bibkey.tJO = 'BaskGvoz2021';
comment.tJO = 'Individual data';

%% set weights for all real data
weights = setweights(data, []);
weights.tL_f = 5 * weights.tL_f;
weights.tL_m = 5 * weights.tL_m;
weights.tWw = 5 * weights.tWw;

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
set2 = {'tWw'}; subtitle2 = {'Wet weight'};
set3 = {'tJO'}; subtitle3 = {'O2 consumption'};
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
