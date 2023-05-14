L = 56;
res =100;

minR = log10(2*L/(pi));
maxR = log10(200);
Rc = [logspace(minR,maxR,res/2)]';

z = Rc.*sin(L./Rc);

%[est_Rc, gof] = fit(z,Rc,'exp1','normalize','on')
[est_Rc, gof] = fit(z,Rc,'fourier5')

figure(1)
plot(est_Rc,z,Rc)
ylabel('Rc')
xlabel('z')