function [PF z Z1 neZ zUpperLim xMaximum f_reflex XLinear ZLinear ZRefLinear NLX NLZ breakpoint alpha_nVector] = IonosphereRayTracing(ParameterArray);
    
    %Physical constants
    K_b=1.38065e-23;%Botlzmann contstant
    m0=1.672622e-27; %Proton mass
    g0=9.8067;% g acceleration at sea level
    q=1.6e-19;%Electron charge
    E0=8.854e-12;%Electric permittivity in vacuum
    m_e=9.10938e-31;%Electron mass
    m=16*2*m0; %O2 mass
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Get parameters from param().m
    % Some fixed values for debug
    %alpha=2*pi/6; %Ray angle in radians (over the vertical)
    %a_i=9.3e-23;% Ionization coefficient
    %a_r=3e-14;% Recombination coefficient
    %a_a=9.3e-23;% Absorttion coefficient
    %T=500 %Temperature in kelvin
    %f=5e6 %Frequency in Hz
    %n0=2.7e25; % Neutral species density at sea level
    %I0=1.5e14;% Radiation intensity at the upper layer of the atmosphere
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    deltaZ=100; %deltaZ in meters
    zUpperLim=6e5; %Maximum z computed in m
    MaxComputedAngle= 89; %Maximum computed angle for step 4 in degrees
    
    %ParameterArrayValues
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    a_a=ParameterArray(1);% Absorttion coefficient
    a_i=ParameterArray(2);% Ionization coefficient
    a_r=ParameterArray(3);% Recombination coefficient
    I0=ParameterArray(4);% Radiation intensity at the upper layer of the atmosphere
    n0=ParameterArray(5); % Neutral species density at sea level
    T=ParameterArray(6);%Temperature in kelvin
    f=ParameterArray(7);%Frequency in Hz
    alpha=ParameterArray(8);%Ray angle in degrees (over the vertical)
    H=K_b*T/(m*g0); %Atmospheric scale height
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    z=0:1000:zUpperLim; %Z testing vector
    ne=@(z_)((a_i/a_r)*n0.*exp(-(z_)/H)*I0.*exp(-H*a_a*n0*exp(-(z_)/H))).^0.5; %Electron density at chapman layer mathematical function
    neZ=ne(z); %Electron density for z vector
    f_p=@(n_)((n_.*((q^2)/(E0*m_e))).^0.5)/(2*pi); %plasma frequency for a given electron density
    f_reflex=f*cos(alpha); %Theoretical frequency of reflexion
    PF=f_p(ne(z));%Plasma frequency for each Z
    Z1=0; %Resolution Value change
    
    %Getting the relexion height and critical height Z1
    for k=0:deltaZ:zUpperLim
        fpi=f_p(ne(k));
        if fpi>(f_reflex/10) && Z1==0
            Z1=k-deltaZ;
        end
        if fpi>f_reflex
            zreflex=k-deltaZ;
            break
        end
    end
    
    Zscaled=[0:floor(Z1/10):Z1,Z1:deltaZ:zreflex]; %Z vector with only 10 samples until Z1 point (linear part), then samples at deltaZ resolution defined previously
    [xz,nzreflexion]=size(Zscaled); %Get the Zscaled vector lenght "nzreflexion" used for iteration later
    xMaximum=2*Zscaled(nzreflexion)*tan(alpha); %Maximum travelled distance (linearly), used to flip the linear ray path    
    
    %Step 3 define the simplified path with origin, reflexion and destination point
    XLinear=[0,zreflex*tan(alpha),2*zreflex*tan(alpha)];
    ZLinear=[0,zreflex,0];
    ZRefLinear=zreflex;
    
    
    %Add a deltaZ lenght with 10 times the resolution to approach better to the critical angle
    ZscaledExtra=[Zscaled,(Zscaled(nzreflexion)+0.1*deltaZ):(0.1*deltaZ):(Zscaled(nzreflexion)+deltaZ)]; 
    [aa,bb]=size(ZscaledExtra);
    x1=zeros(size(ZscaledExtra)); %X coordinate calculate with vatiable refracting index
    Dx1=zeros(size(ZscaledExtra)); %Delta x1 used to flip the x1 vector
    alpha_n=alpha; %initial value for iterative aplha
    alpha_nVector=zeros(1,nzreflexion); %Vector to save angle value (in degrees) for each iteration step (debug)
    nr=@(z) (1-((f_p(ne(z))^2))/(f^2))^0.5; %Refractive index function
    breakpoint=0; %index (iteration step) when alpha_n reaches its criticall value "MaxComputedAngle"
    
    %Loop for each index of ZscaledExtra (1:1:bb)
    for k=1:1:bb
        if k==1
            %For the first index as we don't have a previous value 
            alpha_n=alpha;
            x1(k)=ZscaledExtra(k)*tan(alpha_n);
            Dx1(k)=x1(k);
        else
            %Compute the new incidence angle, if is greater than 90 stop the loop and mark the previous index as the last for the final vector
            alpha_n=asin((nr(ZscaledExtra(k-1))*sin(alpha_n))/(nr(ZscaledExtra(k))));
            if alpha_n>(MaxComputedAngle*pi/180)
                breakpoint=k-1;
                break
            end
            %The difference between the current X and the one previously computed is saved to flip the vector later
            Dx1(k)=(ZscaledExtra(k)-ZscaledExtra(k-1))*tan(alpha_n);
            %New position is computed as the difference between the current and previous Z times the tangent of alpha for its corresponding difraction index
            x1(k)=x1(k-1)+Dx1(k);            
        end
        alpha_nVector(k)=180*alpha_n/pi; %Vector with all the diffraction angles computed for each step
    end
    %Trim the final vectors to the breakpoint and then flip Z and X
    Dx11=Dx1(1:breakpoint);
    flipDX=flip(Dx11);
    x11=x1(1:breakpoint);
    x2=zeros(1,breakpoint);
    alpha_nVector=(alpha_nVector(1:breakpoint));
    for k=1:1:breakpoint
        if k==1
            x2(k)=x11(breakpoint)+flipDX(k);
        else
            x2(k)=x2(k-1)+flipDX(k-1);
        end
    end
    
    XCombinedNL=[x11,x2]; %X combined non linear
    ZCombinedNL=[ZscaledExtra(1:breakpoint),flip(ZscaledExtra(1:breakpoint))];%Z combined non linear
    NLX=XCombinedNL;
    NLZ=ZCombinedNL;    
   
end