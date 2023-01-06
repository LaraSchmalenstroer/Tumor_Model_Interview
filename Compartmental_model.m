load("parameters.mat")
load('initial_conditions.mat')

m = sbiomodel('Compartmental model');
c1 = addcompartment(m, 'tumor');
c2 = addcompartment(m, 'central');

% add species and their initial values to the model 
s1 = addspecies(c1, 'T1', 'InitialAmount', T1);
s2 = addspecies(c1, 'T2', 'InitialAmount', T2);
s3 = addspecies(c1, 'E1', E1);
s4 = addspecies(c1, 'E2', E2);
s5 = addspecies(c2, 'L', 5e7);

% add parameters and their values
p1 = addparameter(m, 'g1', g1);
p2 = addparameter(m, 'g2', g2);
p3 = addparameter(m, 'K1', K1);
p4 = addparameter(m, 'K2', K2);
p5 = addparameter(m, 'c12', c12);
p6 = addparameter(m, 'c21', c21);
p7 = addparameter(m, 'a11', a11);
p8 = addparameter(m, 'a12', a12);
p9 = addparameter(m, 'a21', a21);
p10 = addparameter(m, 'p1', p_1);
p11 = addparameter(m, 'd1', d1);
p12 = addparameter(m, 'd2', d2);
p13 = addparameter(m, 'e1', e1);
p14 = addparameter(m, 'e2', e2);
p15 = addparameter(m, 'r1', r1);
p16 = addparameter(m, 'r2', r2);
p17 = addparameter(m, 'r3', r3);
p18 = addparameter(m, 's1', s_1);
p19 = addparameter(m, 's2', s_2);

% add equations
r1 = addrule(m, 'T1=g1*T1*(1-(T1/K1))-a11*E1*T1-a12*E2*T1-c12*T1*T2', ...
    'RuleType', 'rate');
r2 = addrule(m, 'T2=g2*T2*(1-(T2/K2))-a21*E1*T2-c21*T1*T2', ...
    'RuleType', 'rate');
r3 = addrule(m, 'E1=p1-d1*E1-e1*(T1+T2)*E1+(r1*(T1+T2))/(s1+T1+T2)*E1', ...
    'RuleType', 'rate');
r4 = addrule(m, 'E2=-d2*E2-e2*T1*E2+(r2*T1)/(s2+T1)*E2+r3*E1*(T1+T2)', ...
    'RuleType', 'rate');
r5 = addrule(m, 'L=e1*T1+K1*E1-K2*E1', 'RuleType', 'rate');

csObj = m.addconfigset('newStopTimeConfigSet');
csObj.StopTime = 300;
[t,sd,species] = sbiosimulate(m, csObj);

%plot species dynamics vs. time
plot(t, sd);
legend(species)
xlabel('Time');
ylabel('Amount species');

%plot tumor vs. E1 cells
% tumor_total = sd(:,1) + sd(:,2);
% plot(sd(:,3), tumor_total);
% legend('tumor total')
% xlabel('Amount E1 cells');
% ylabel('Amount tumor cells');
