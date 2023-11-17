%[PF Z0 neZ zUpperLim xMaximum f_reflex XLinear ZLinear ZRefLinear NLX NLZ breakpoint]=
[aa ai ar I0 n0 T f alpha]=get_parameters(980131);
[PF Z0 Z1 neZ zUpperLim xMaximum f_reflex XLinear ZLinear ZRefLinear NLX NLZ breakpoint alpha_nVector] = IonosphereRayTracing([aa ai ar I0 n0 T f ((alpha*pi)/180)]);

%Refraction angle vs altitude FIG8
fig8=figure(8);
plot(alpha_nVector,NLZ(1:breakpoint)/1e3,"-*")
yline(Z1/1e3,"-.",["z(f_r_e_f_l_x/10): " + num2str(Z1/1e3,4)+ "km"])
xlabel('alpha [deg]')
ylabel('z [km]')
title("Refraction angle vs altitude")

fig1=figure(1);
semilogx(neZ,Z0/1000) %SemiLog in x plot of electron density vs altitude (z in km) FIGURE1
xlabel('ne [m^-^3]')
ylabel('z [km]')
title("Electron density vs altitude")
fig2=figure(2);
semilogx(PF/1000,Z0/1000); %SemiLog in x plot of flasma frequency (kHz) vs altitude (z in km), 
% Vertical line of emission frequency FIGURE2
xlim([0.1,2e8/1e3])
ylim([0,zUpperLim/1e3])
xline(f_reflex/1e3,"-.",["Reflexion frequency: " num2str(f_reflex/1e3,4)+ "kHz"],LabelOrientation="horizontal") %Plot vertical line at reflexion frequency
yline(ZRefLinear/1e3,"-.",["Reflexion altitude: " + num2str(ZRefLinear/1e3,4)+ "km"],"LabelHorizontalAlignment","left") %Plot horizontal line at reflexion altitude
xlabel('f_p [kHz]')
ylabel('z [km]')
title("Plasma frequency vs altitude")

fig3=figure(3);
plot(XLinear/1e3,ZLinear/1e3)
%Plot of the coordinates of a linear ray path (both axis in km) FIGURE3
yline(ZRefLinear/1e3,"-.",["Reflexion altitude: " + num2str(ZRefLinear/1e3,4)+ "km"]) %Plot horizontal line at reflexion altitude
xlabel('x [km]')
ylabel('z [km]')
xlim([0.1,(xMaximum/1e3)+50]) %End the X axis at xMaximum value+50km
title("Linear aproximation ray path")

fig4=figure(4);
%the reflection shall be computed with the difference of the flipped vector Zi
plot(NLX/1e3,NLZ/1e3)
yline(NLZ(breakpoint)/1e3,"-.",["Reflexion altitude: "+num2str(NLZ(breakpoint)/1e3,4)+"km"]) %Plot horizontal line at reflexion altitude
xlabel('x [km]')
ylabel('z [km]')
xlim([0,(NLX(2*breakpoint)/1e3)+50]) %End the X axis at xMaximum value + 50km
title("Variable Refracive index aproximation ray path")

fig5=figure(5);
%Comparisson Figure
plot(XLinear/1e3,ZLinear/1e3,NLX/1e3,NLZ/1e3)
%Plot of the coordinates of a linear ray path (both axis in km) FIGURE3
yline(ZRefLinear/1e3,"-.",["Reflexion altitude (linear): "+num2str(ZRefLinear/1e3,4)+ "km"],LabelVerticalAlignment="bottom") %Plot horizontal line at reflexion altitude (linear)
yline(NLZ(breakpoint)/1e3,"-.",["Reflexion altitude: "+num2str(NLZ(breakpoint)/1e3,4)+ "km"]) %Plot horizontal line at reflexion altitude
xlabel('x [km]')
ylabel('z [km]')
xlim([0,(NLX(2*breakpoint)/1e3)+50]) %End the X axis at xMaximum value +50km
title("Linear vs Variable Refracive index aproximation")

%STEP 5
[PF1 Z01 Z1 neZ1 zUpperLim1 xMaximum1 f_reflex1 XLinear1 ZLinear1 ZRefLinear1 NLX1 NLZ1 breakpoint1] = IonosphereRayTracing([aa ai ar I0 n0 T+100 f ((alpha*pi)/180)]);

fig6=figure(6);
%Plot of fig4 with original parameters and with T1=T+100K
plot(NLX/1e3,NLZ/1e3,NLX1/1e3,NLZ1/1e3)
yline(NLZ(breakpoint)/1e3,"-.",["Reflexion altitude: "+num2str(NLZ(breakpoint)/1e3,4)+"km"]) %Plot horizontal line at reflexion altitude
yline(NLZ1(breakpoint1)/1e3,"-.",["Reflexion altitude +100K: "+ num2str(NLZ1(breakpoint1)/1e3,4)+ "km"])
xlabel('x [km]')
ylabel('z [km]')
title("Variable Refracive index aproximation: "+num2str(T,4)+"K vs "+num2str(T+100,4)+"K")

fig9=figure(9);
semilogx(PF/1000,Z0/1000,PF1/1000,Z0/1000); %SemiLog in x plot of flasma frequency (kHz) vs altitude (z in km), 
% Vertical line of emission frequency FIGURE2
xlim([0.1,2e8/1e3])
ylim([0,zUpperLim/1e3])
yline(ZRefLinear1/1e3,"-.",["Reflexion altitude +100K: " + num2str(ZRefLinear1/1e3,4)+ "km"],"LabelHorizontalAlignment","left") %Plot horizontal line at reflexion altitude
xline(f_reflex/1e3,"-.",["Reflexion frequency: " num2str(f_reflex/1e3,4)+ "kHz"],LabelOrientation="horizontal") %Plot vertical line at reflexion frequency
yline(ZRefLinear/1e3,"-.",["Reflexion altitude: " + num2str(ZRefLinear/1e3,4)+ "km"],"LabelHorizontalAlignment","left",LabelVerticalAlignment="bottom") %Plot horizontal line at reflexion altitude
xlabel('f_p [kHz]')
ylabel('z [km]')
title("Plasma frequency vs altitude")

%Step 7
[PF2 Z02 Z1 neZ2 zUpperLim2 xMaximum2 f_reflex2 XLinear2 ZLinear2 ZRefLinear2 NLX2 NLZ2 breakpoint2] = IonosphereRayTracing([aa ai ar I0*10 n0 T+100 f ((alpha*pi)/180)]);

fig7=figure(7);
%Plot of fig4 with original parameters and with 10 times I0
plot(NLX/1e3,NLZ/1e3,NLX2/1e3,NLZ2/1e3)
yline(NLZ(breakpoint)/1e3,"-.",["Reflexion altitude: "+num2str(NLZ(breakpoint)/1e3,4)+ "km"]) %Plot horizontal line at reflexion altitude
yline(NLZ2(breakpoint2)/1e3,"-.",["Reflexion altitude 10*I_0: "+num2str(NLZ2(breakpoint2)/1e3,4)+ "km"],LabelVerticalAlignment="bottom")
xlabel('x [km]')
ylabel('z [km]')
xlim([0,(NLX(2*breakpoint)/1e3)+50]) %End the X axis at xMaximum value +50km
title("Variable Refracive index aproximation I_0: "+num2str(I0,'%.2e')+" vs "+num2str(10*I0,'%.2e'))


fig10=figure(10);
semilogx(PF/1000,Z0/1000,PF2/1000,Z0/1000); %SemiLog in x plot of flasma frequency (kHz) vs altitude (z in km), 
% Vertical line of emission frequency FIGURE2
xlim([0.1,2e8/1e3])
ylim([0,zUpperLim/1e3])
yline(ZRefLinear2/1e3,"-.",["Reflexion altitude 10*I0: " + num2str(ZRefLinear2/1e3,4)+ "km"],"LabelHorizontalAlignment","left",LabelVerticalAlignment="bottom") %Plot horizontal line at reflexion altitude
xline(f_reflex/1e3,"-.",["Reflexion frequency: " num2str(f_reflex/1e3,4)+ "kHz"],LabelOrientation="horizontal") %Plot vertical line at reflexion frequency
yline(ZRefLinear/1e3,"-.",["Reflexion altitude: " + num2str(ZRefLinear/1e3,4)+ "km"],"LabelHorizontalAlignment","left") %Plot horizontal line at reflexion altitude
xlabel('f_p [kHz]')
ylabel('z [km]')
title("Plasma frequency vs altitude")

str = ['x max: Linear aprox    = ' num2str(XLinear(3)/1e3,3) ' | Realistic aprox =  ' num2str(NLX(2*breakpoint)/1e3,3) ' (km) || z max: Linear aprox    = ' num2str(ZLinear(2)/1e3,4) '  | Realistic aprox =  ' num2str(NLZ(breakpoint)/1e3,4) ' (km)'];
disp(str);
str = ['x max: original values = ' num2str(NLX(2*breakpoint)/1e3,3) ' | T+100 values    =  ' num2str(NLX1(2*breakpoint1)/1e3,3) ' (km) || z max: original values = ' num2str(NLZ(breakpoint)/1e3,4) ' | T+100 values    =  ' num2str(NLZ1(breakpoint1)/1e3,4) '  (km)'];
disp(str);
str = ['x max: original values = ' num2str(NLX(2*breakpoint)/1e3,3) ' | 10*I0 values    =  ' num2str(NLX2(2*breakpoint2)/1e3,3) ' (km) || x max: original values = ' num2str(NLZ(breakpoint)/1e3,4) ' | 10*I0 values    =  ' num2str(NLZ2(breakpoint2)/1e3,4) ' (km)'];
disp(str);