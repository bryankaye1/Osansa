astock = 250;
alabel = 125;
donor = 200;
ex = 24;


deV = (20*ex + 1*donor)/21;

a1V = (10*ex + 1*astock)/11;
a1_labV = alabel/11;
a1_labN= a1_labV/a1V;

a2V = (20*ex + 20*a1V)/42;
a2_labV = a1_labV/(42/20);
a2_lab = a2_labV/a2V;

a3V = (15*ex + 15*a2V)/32;
a3_labV = a2_labV/(32/15);
a3_lab = a3_labV/a3V;

m1V = (a1V+deV)/2; 
m2V = (a2V+deV)/2;
m3V = (a3V+deV)/2;

m1_labV = a1_labV/2;
m1_labN = m1_labV/m1V;

m2_labV = a2_labV/2;
m2_labN = m2_labV/m2V;

m3_labV = a3_labV/2;
m3_labN = m3_labV/m3V;

fprintf('m1_Afrac %1.3f m2_Afrac is %1.3f m3_Afrac is %1.3f \n', 200*m1_labN,200*m2_labN,200*m3_labN);

aest = [m1_labN,m1_labN,m1_labN,  m2_labN,m2_labN,m2_labN, m3_labN,m3_labN,m3_labN];
af = [.304,.304,.293,.186,.213,.197,.168,.14,.184];

figure(2); clf;scatter(aest,af);


aestm = [m1_labN,m2_labN,m3_labN ];
afm = [mean(af(1:3)),mean(af(4:6)),mean(af(7:9))];
figure(3); clf;scatter(aestm,afm);





