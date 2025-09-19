function proj4_v2()
    N = 10000;
    M = 50000;
    r = 1/100;
    mu = 1/10;
    T = 100;
    n = 200;
    h = T/n;
    t = linspace(0, T, n+1);


    x0 = 0.01; 
    y0 = 0.01;
    yinit = [x0; y0];
    t0 = 0;


    aval = [0.2, 0.3, 0.4];
    bval = [0.1, 0.2, 0.3];

    clf;

    for i = 1:length(aval)
        for j = 1:length(bval)
            a = aval(i);
            b = bval(j);


            function dydt = f(t, y)
                x = y(1);
                m = y(2);
                dx = (a*b*M/N)*m*(1-x) - r*x;
                dm = a*x*(1-m) - mu*m;
                dydt = [dx; dm];
            endfunction


            y = ode(yinit, t0, t, f);


            subplot(length(aval), length(bval), (i-1)*length(bval)+j);
            plot(t, y(1,:), 'g-', t, y(2,:), 'b-');
            xlabel("Temps");
            ylabel("Population infectée");
            legend("Humains (x)", "Moustiques (y)", 4);
            title("a="+string(a)+", b="+string(b));
        end
    end
endfunction

proj4_v2();
