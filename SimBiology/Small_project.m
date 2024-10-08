% create model 
m = sbiomodel('small_project');

% add species and their initial values to the model 
s1 = addspecies(m, 'T1', 'InitialAmount', 8e7);
s2 = addspecies(m, 'T2', 'InitialAmount', 2e7);
s3 = addspecies(m, 'E1', 1.1e7);
s4 = addspecies(m, 'E2', 0);

% add parameters and their values
p1 = addparameter(m, 'g1', 5.14e-1);
p2 = addparameter(m, 'g2', 0.35*p1.get('Value'));
p3 = addparameter(m, 'K1', 5e8);
p4 = addparameter(m, 'K2', p3.get('Value'));
p5 = addparameter(m, 'c12', 1.1e-9);
p6 = addparameter(m, 'c21', 1.5*p5.get('Value'));
p7 = addparameter(m, 'a11', 1.1e-7);
p8 = addparameter(m, 'a12', 1.1e-10);
p9 = addparameter(m, 'a21', p7.get('Value'));
p10 = addparameter(m, 'p1', 1.3e4);
p11 = addparameter(m, 'd1', 4.12e-2);
p12 = addparameter(m, 'd2', 2e-2);
p13 = addparameter(m, 'e1', 3.42e-10);
p14 = addparameter(m, 'e2', 3.42e-10);
p15 = addparameter(m, 'r1', 1.24e-1);
p16 = addparameter(m, 'r2', 1.24e-3);
p17 = addparameter(m, 'r3', 1.1e-7);
p18 = addparameter(m, 's1', 2.02e7);
p19 = addparameter(m, 's2', 2.02e7);

% add equations
r1 = addrule(m, 'T1=g1*T1*(1-(T1/K1))-a11*E1*T1-a12*E2*T1-c12*T1*T2', ...
    'RuleType', 'rate');
r2 = addrule(m, 'T2=g2*T2*(1-(T2/K2))-a21*E1*T2-c21*T1*T2', ...
    'RuleType', 'rate');
r3 = addrule(m, 'E1=p1-d1*E1-e1*(T1+T2)*E1+(r1*(T1+T2))/(s1+T1+T2)*E1', ...
    'RuleType', 'rate');
r4 = addrule(m, 'E2=-d2*E2-e2*T1*E2+(r2*T1)/(s2+T1)*E2+r3*E1*(T1+T2)', ...
    'RuleType', 'rate');

csObj = m.addconfigset('newStopTimeConfigSet');
csObj.StopTime = 300;
[t,sd,species] = sbiosimulate(m, csObj);

%plot species dynamics vs. time
plot(t, sd);
legend(species)
xlabel('Time');
ylabel('Amount species');

%plot tumor vs. E1 cells
tumor_total = sd(:,1) + sd(:,2);
plot(sd(:,3), tumor_total);
legend('tumor total')
xlabel('Amount E1 cells');s
ylabel('Amount tumor cells');
