function [pa_c,daa_c] = pure_pa_diffusion_coef(wave,enviro,energy_tick)
    
    %It's to compute the pure pitch angle diffusion coef averaged over one
    %bounce orbit as suggested by Lawrence Lyons in the 1970s.
    %This algorithm doesn't assume the wave frequency distribution to Gaussian
    %or any specific form.It feeds on the measured wave power spectral.
    
    % it's to compute pure pitch angle diffusion, more specifically, low
    % frequency whistler mode such as Hiss interaction with radiation belt
    % electrons. Laudau resonance at n = 0 is ignored given that
    % diffuison at small pitch angles inside losscone is interested. Landau
    % resonance diffusion coef is several orders of magnitude smaller than
    % cycltron resonance diffusion coef inside losscone.
    
    % The unit system is CGS.
    
    % by Leo Zheng , Oct 27th, 2017.
    
    %inputs: wave, wave parameters
    %           .f is wave frequency in Hz, frequency stamp of .bw2
    %           .bw2 is a vector of wave spectral density in Gauss2/Hz
    %           .Xm is tan(central normal angle)
    %           .Xw is the width of tan(normal angle)
    %           .wlc, lower end cut off frequency in rad/s 
    %       enviro.L is L shell, enviro.mlt is magnetic locall time,
    %       enviro.density is cold plasma density in #/cm^3 
    %       energy_tick, in keV, default values are from
    %        CRRES MEA energy bins
    
    
    %outputs:daa_c, bounce averaged cyclotron resonance diffusion coefs,
    %          unit 1/s, a 2-d matrix with dimension of len(energy_tick) by len(pa_c)
    %         pa_c, pitch angle stamps for daa_c in ascending order
    %          division of pitch angles within a loss cone, in degree
    
    %calculate size of loss cone in degree
    loss_cone = cal_LossCone(enviro.L);
    pa_c = linspace(0,loss_cone,6);
    energy_tick = sort(energy_tick);
    nE = length(energy_tick);
    nPa = length(pa_c);
    daa_c = zeros(nE,nPa);
    nmin = -5; % highest negative harmonic
    nmax = 5; % hightest positive harmonic

    for i = 1: nE
        for j = 1:nPa
            particle.pa = pa_c(j); % degree
            particle.kE = energy_tick(i); % keV 

            disp('iteration:')
            disp([i,j])
            
           % compute bounce-averaged pure pitch angle diffusion coefs at
           % given particle energy and pitch angle
            daa_c(i,j) = compute_daa_c(particle, enviro, wave,nmin,nmax);
            
        end
        
    end
    

end

function daa = compute_daa_c(particle, enviro, wave,nmin,nmax)
    % compute bounce averaged pure pitch angle diffusion coefs
    
    % fundamental constants in CGS unit
    e = 4.8032e-10; % elementary charge in esu
    me = 9.10938356e-28 ; % electron rest mass in gram
    c = 2.998e10; % speed of light in cm/sec
    eqB0 = 0.31; % equatorial magnetic field on the ground in Gauss
    
    %initialize parameters 
    % particle and environment
    para.pa = particle.pa * pi/180; % convert degree to radian
    para.kE = particle.kE;
    para.L = enviro.L;
    % wave
    para.f = wave.f;
    para.bw2 = wave.bw2;
    para.Xw = wave.Xw;
    para.Xm = wave.Xm;
    para.wlc = wave.wlc;
    % initialize depdent parameters    
    para.m2c2 = (me*c)^2;
    para.gamma = 1 + particle.kE/511; % relativistic factor
    para.eqB =  eqB0/(enviro.L)^3; %equatorial magnetic field
    para.eqGyoFreq = e * para.eqB/(me * c); % equatorial gyrofrequency
    para.wpe2 = 4*pi*enviro.density*e^2/me; % plasma frequency square
    para.p = (sqrt((1 + particle.kE/511)^2 - 1))*me*c; % momentum of particle
    
    para.bounceIntegral = 1.3 - 0.56*sin(para.pa); % pitch angle function for bounce period
    
    % compute mirroring latitude
    para.mirrorLat = get_mirrorLat(para);
    
    daa = pi*e^2*(para.eqGyoFreq)^2/8/para.gamma/(para.p)^2/para.wpe2;
    daa = daa/cos(para.pa)^2/para.bounceIntegral;
    
    lambda_max = para.mirrorLat;
    I = integral(@(x)LatFun4DaaCyc(para,x,nmin,nmax),0,lambda_max);
    daa = daa*I;
    
    
end

% compute mirroring latitude 
function mirrorLat = get_mirrorLat(para)   
            % h(lat) = 1/sin(pa)^2; h(lat) = (1+3*sin(lat)^2)^0.5/cos(lat)^6
            % x = sin(lat)^2, x [0,1]
            % tranformed to sin(pa)^(2/3) * (1+3*x)^(1/6) = 1-x, so it can
            % handle pa = 0 case and returns only one root rather than 6
            % identical roots.
            k = sin(para.pa).^(2./3.);
            syms x 
            s = vpasolve(1-x == k*(1+3*x)^(1./6),x,[0,1]);
            s2 = 1-1./para.L; % invariant latitude cos(invarLat)^2 = 1/L, x = sin(lat)^2, the dipole field latitude cutoff is invariant latitude rather than 90degree 
            s = min(double(s),s2);
            mirrorLat = asin(sqrt(s)); % in radian   
            
end

function fun = LatFun4DaaCyc(para,lat,nmin,nmax)
            h_lambda = ((1+3*sin(lat).^2).^(0.5))./(cos(lat)).^6;
            fun = (cos(lat).^7).*(h_lambda.^2)./I_normal_angle(para)./sqrt(1-h_lambda*(sin(para.pa)^2));
            fun = fun.*I_harmonics(para,nmin,nmax,lat);
end
function I = I_normal_angle(para)
    % this integral is twice of the normal angle integral I in Appendix C
    % of Lyons,Thorne and Kennel 1971. The same I is part of W1 in Lyons
    % 1972.
    I = integral(@(x)fun4I_normal_angle(para,x),0,inf);
    
end
function f = fun4I_normal_angle(para,x)
    f = x./((1+x.^2).^(0.75)); f = f.*wave_normal_anlge_model(para,x);
end

function model = wave_normal_anlge_model(para,x)
    arg = -(x - para.Xm).^2./(para.Xw)^2;
    model =  exp(arg);
end

function I = I_harmonics(para,nmin,nmax,lat)
   if(nmin > nmax) 
        disp('n range is wrong. nmin should not be larger than nmax')
        return;
   end
   if(nmin * nmax == 0)
        disp('no zero order harmonic for cyclotron resonance. Input non-zero n please.')
        return;
   end
   I = 0;
   for n = nmin:nmax
       if(n ~= 0)
          I = I + abs(n)*I_n(para,n,lat); 
       end
   end
end

function I = I_n(para,n,lat)
    h_lambda = ((1+3*sin(lat).^2).^(0.5))./(cos(lat)).^6;
    % compute lower limit of x(tan(theta)) given lower cut off frequency
    % and k-w dispersion relation
    dispersion = para.wlc*para.wpe2*(para.p)^2*(1 - h_lambda *(sin(para.pa))^2);
    dispersion = dispersion/n^2/para.m2c2/(para.eqGyoFreq)^3./(h_lambda).^3;
    dispersion = dispersion.^2;
    xmin = zeros(size(dispersion));
    idx = dispersion > 1;
    xmin(idx) = sqrt(dispersion(idx) - 1);
    
    I = arrayfun(@(lat,xmin)integral(@(x)fun4I_n(para,n,lat,x),xmin,inf,'ArrayValued',true),lat,xmin);
    
end
function f = fun4I_n(para,n,lat,x)
   f = x./((1+x.^2).^0.25);
   f = f.*wave_power_spec(para,n,lat,x);
   f = f.*wave_normal_anlge_model(para,x);
   f = f.*Bessel_n(para,n,lat,x);
end

% given data of wave power spectra on frequency domain, find the power
% spectra at frequencies corresponding to tan(theta) and resonance k||, theta is wave normal
% angle, k|| is the parallel component of wave vector at resonance. Greatly
% simplified cold plamsma Hiss dispersion relation is used (A6 in Lyons,Thorne and Kennel 1971)
function bw2 = wave_power_spec(para,n,lat,x)
    freq = para.f; % in Hz
    bw2 = para.bw2;  % in Gauss2/Hz
    bw2(isnan(bw2)) = 0;
    % compute frequency given x and k-w disperion relation, x is tan(theta)
    h_lambda = ((1+3*sin(lat).^2).^(0.5))./(cos(lat)).^6;
    freq_x = n^2*para.m2c2*((para.eqGyoFreq)^3)/para.wpe2/(para.p)^2;
    freq_x = freq_x*h_lambda.^3./(1 - sin(para.pa)^2 * h_lambda);
    freq_x = freq_x*sqrt(1+x.^2);
    freq_x = freq_x/2/pi; % convert angular frequency to frequency
    
    if(freq_x < min(freq) || freq_x > max(freq))
        bw2_int = 0;
    else
        bw2_int = interp1(freq,bw2,freq_x);
    end
    
    bw2 = bw2_int;
    
end

function bes = Bessel_n(para,n,lat,x)
    % x = tan(theta)
    coef_pos = 1 + (1+x.^2).^(-0.5);
    coef_neg = 1 - (1+x.^2).^(-0.5);
    a0 = para.pa;
    %arg = sin(a0)^2 * Daa.hLamda(lamda);
    h_lamdda = ((1+3*sin(lat).^2).^(0.5))./(cos(lat)).^6;
    arg = sin(a0)^2 * h_lamdda;
    arg = sqrt(arg/(1-arg));
    arg = -n*arg*x ;
    bes = (coef_pos.*besselj(n+1,arg) + coef_neg.*besselj(n-1,arg)).^2;
    
end



