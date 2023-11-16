function [PF Z0 neZ zUpperLim xMaximum f_reflex XLinear ZLinear ZRefLinear NLX NLZ breakpoint] = IonosphereRayTracing(ParameterArray)
    
    K_b=1.38065e-23;%Botlzmann contstant
    m0=1.672622e-27; %Proton mass
    g0=9.8067;% g acceleration at sea level
    Ro=1.225;%atmospheric density at sea level
    q=1.6e-19;%Electron charge
    E0=8.854e-12;
    m_e=9.10938e-31;
    zUpperLim=6e5 %Maximum z computed in m
    %Get parameters from param().m
    % 31/01/1998
    alpha=2*pi/6; %Ray angle in radians (over the vertical)
    a_i=9.3e-23;% Ionization coefficient
    a_r=3e-14;% Recombination coefficient
    a_a=9.3e-23;% Absorttion coefficient
    T=500 %Temperature in kelvin
    f=5e6 %Frequency in Hz
    deltaZ=100 %deltaZ in meters
    n0=2.7e25; % Neutral species density at sea level
    I0=1.5e14;% Radiation intensity at the upper layer of the atmosphere
    m=16*2*m0; %O2 mass
    MaxComputedAngle= 89 %Maximum computed angle for step 4 in degrees
    
    %ParameterArrayValues
    a_a=ParameterArray(1)
    a_i=ParameterArray(2)
    a_r=ParameterArray(3)
    I0=ParameterArray(4)
    n0=ParameterArray(5)
    T=ParameterArray(6)
    f=ParameterArray(7)
    alpha=ParameterArray(8)
    H=K_b*T/(m*g0); %Atmospheric scale height
    z=0:deltaZ:zUpperLim; %Z testing vector
    Z0=z
    ne=@(z_)((a_i/a_r)*n0.*exp(-(z_)/H)*I0.*exp(-H*a_a*n0*exp(-(z_)/H))).^0.5; %Electron density at chapman layer mathematical function
    neZ=ne(z)
    f_p=@(n_)((n_.*((q^2)/(E0*m_e))).^0.5)/(2*pi); %plasma frequency for a given electron density
    f_reflex=f*cos(alpha) %Theoretical frequency of reflexion
    PF=f_p(ne(z));
    [dim, steps]=size(z); %CLEAN
    Z1=0; %Resolution Value change
    
    %Try to get electron density critical values by math instead of bruteforce
    for k=0:deltaZ:zUpperLim
        fpi=f_p(ne(k));
        if fpi>(f/10) && Z1==0
            Z1=k-deltaZ-1
        end
        if fpi>f_reflex
            zreflex=k-deltaZ-1
            break
        end
    end
    

    
    Zscaled=[0:floor(Z1/5):Z1,Z1:deltaZ:zreflex]; %Z vector with only 5 samples until Z1 point (linear part), then samples at deltaZ resolution defined previously
    [xz,nzreflexion]=size(Zscaled); %Get the Zscaled vector lenght "nzreflexion" used for iteration later
    Xscaled=zeros(1,(2*nzreflexion)); %Define a X vector with twice the size of ZScaled
    xMaximum=2*Zscaled(nzreflexion)*tan(alpha); %Maximum travelled distance (linearly), used to flip the linear ray path
    for k=1:1:nzreflexion
        %Compute the x value for each Z sampled (as x is mirrored is x(k) and xMaximun-x(k) are computed at the same time)
        Xlinear=Zscaled(k)*tan(alpha);
        Xscaled(k)=Xlinear;
        Xscaled((2*nzreflexion)-k+1)=xMaximum-Xlinear;
    end
    ZCombined=[Zscaled,flip(Zscaled)]; %
    XLinear=Xscaled;
    ZLinear=ZCombined;
    ZRefLinear=zreflex;
    
    
    Extra=(Zscaled(nzreflexion)+0.1*deltaZ):(0.1*deltaZ):(Zscaled(nzreflexion)+deltaZ)
    ZscaledExtra=[Zscaled,Extra] %Add a deltaZ lenght with 10 times the resolution to approach better to the critical angle
    [aa,bb]=size(ZscaledExtra)
    x1=zeros(size(ZscaledExtra)); %X coordinate calculate with vatiable refracting index
    Dx1=zeros(size(ZscaledExtra)); %Delta x1 used to flip the x1 vector
    alpha_n=alpha; %initial value for iterative aplha
    alpha_nVector=zeros(1,nzreflexion); %Vector to save angle value (in degrees) for each iteration step
    nr=@(z) (1-((f_p(ne(z))^2))/(f^2))^0.5; %Refractive index function
    breakpoint=0 %index (iteration step) when alpha_n reaches its criticall value "MaxComputedAngle"
    
    for k=1:1:bb
        if k==1
            alpha_n=alpha
            x1(k)=ZscaledExtra(k)*tan(alpha_n);
            Dx1(k)=x1(k);
        else
            nr(ZscaledExtra(k));
            (nr(ZscaledExtra(k-1)));
            alpha_n=asin((nr(ZscaledExtra(k-1))*sin(alpha_n))/(nr(ZscaledExtra(k))));
            if alpha_n>(MaxComputedAngle*pi/180)
                breakpoint=k-1
                break
            end
            x1(k)=x1(k-1)+(ZscaledExtra(k)-ZscaledExtra(k-1))*tan(alpha_n);
            Dx1(k)=x1(k)-x1(k-1);
            
        end
        alpha_nVector(k)=180*alpha_n/pi; 
    end
    Dx11=Dx1(1:breakpoint);
    flipDX=flip(Dx11);
    x11=x1(1:breakpoint);
    x2=zeros(1,breakpoint);
    %alpha_nVector=(alpha_nVector(1:breakpoint));
    for k=1:1:breakpoint
        if k==1
            x2(k)=x11(breakpoint)+flipDX(k);
        else
            x2(k)=x2(k-1)+flipDX(k-1);
        end
    end
    
    XCombinedNL=[x11,x2]; %X combined non linear
    ZCombinedNL=[ZscaledExtra(1:breakpoint),flip(ZscaledExtra(1:breakpoint))];%Z combined non linear
    NLX=XCombinedNL
    NLZ=ZCombinedNL

   
end