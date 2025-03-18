clc;clear;close all;

mat_names = ["Fe 430", "Fe 590", "C 45 N", "C 60 N",...
    "34 Cr 4 V", "42 CrMo 4 V", "16 MnCr 5", "15 CrNi 6"];

% module of elasticity and material standards (no price):
% https://www.michael-smith-engineers.co.uk/mse/uploads/resources/useful-info/General-Info/MATERIAL-GRADE-COMPARISON-TABLE-for-Web.pdf
% material (diameter) (grade) (related link)
% Fe 430: (726mm) (S275) https://matweb.com/search/DataSheet.aspx?MatGUID=c2ba59bb365942a7b6da46f1cee370b8  
% Fe 590: (629mm) (S355) https://matweb.com/search/DataSheet.aspx?MatGUID=1dc0414bd1ea4061a5dc09382c455e2a
% C 45 N: (622mm) (AISI/ASTM 1045) https://matweb.com/search/DataSheet.aspx?MatGUID=2ca9b42e83894e8a8a61385fd7da63ae
% C 60 N: (579mm) (SAE/ASTM 1060) https://matweb.com/search/DataSheet.aspx?MatGUID=0a471605c1324daa910855e54a21fab3
% 34 Cr 4 V: (510mm) (AISI 5132/ 1.7033) https://matweb.com/search/DataSheet.aspx?MatGUID=4877d405464f448a96786c8cbd00d3b5
% 42 CrMo 4 V: (487mm) (ASTM A322) https://matweb.com/search/DataSheet.aspx?MatGUID=38108bfd64c44b4c9c6a02af78d5b6c6
% 16 MnCr 5: (285mm) (AISI 5115 / 1.7131) https://matweb.com/search/DataSheet.aspx?MatGUID=2ab813ffa05d40329dffe0ee7f58b5de
% 15 CrNi 6: (277mm) https://matweb.com/search/DataSheet.aspx?MatGUID=9ab3bf332758468ab36010790bd94349 ?

% alibaba pricing:
% Fe 430: https://www.alibaba.com/product-detail/S355JR-ASTM-1045-1055-Cold-Drawn_1601262208795.html?spm=a2700.galleryofferlist.p_offer.d_title.5c4b13a0X6sAvT&s=p
% 16 MnCr 5: https://www.alibaba.com/product-detail/15-Crmn-Smnc420-18-X-16_1601026887161.html?spm=a2700.galleryofferlist.normal_offer.d_title.43d313a0pLzJnI

price_per_kg = [1.5 1.5 1.8 1.8 4 3.5 4.5 5]; % [dollar/kg] very rough estimates from chatGPT, based on rounded up round bar
density = [7.85 7.85 7.80 7.80 7.85 7.85 7.85 7.80]*1e3; % [g/cc] -> [kg/m^3]
prices_m3 = (price_per_kg) .* density;
material_sum_list = [8.5329e+04 5.5350e+04 5.3616e+04 4.2779e+04 2.9276e+04 2.5489e+04 5.1830e+03 4.8950e+03]*1e-6; % [cm^3] -> [m^3]
price_list = material_sum_list .*  prices_m3;
weight = material_sum_list .* density;

% Create the table
material_table = table(mat_names', price_list', weight', 'VariableNames', {'Material Name', 'Price estimate ($)', 'weight (kg)'});
disp(material_table);