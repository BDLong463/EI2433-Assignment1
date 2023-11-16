%[PF Z0 neZ zUpperLim xMaximum f_reflex XLinear ZLinear ZRefLinear NLX NLZ breakpoint]=
[aa ai ar I0 n0 T f alpha]=get_parameters(980131)
[PF Z0 neZ zUpperLim xMaximum f_reflex XLinear ZLinear ZRefLinear NLX NLZ breakpoint] = IonosphereRayTracing([aa ai ar I0 n0 T f ((alpha*pi)/180)])

fig1=figure(1);

semilogx(neZ,Z0/1000) %SemiLog in x plot of electron density vs altitude (z in km) FIGURE1

fig2=figure(2);
semilogx(PF/1000,Z0/1000) %SemiLog in x plot of flasma frequency (kHz) vs altitude (z in km), 
% Vertical line of emission frequency FIGURE2
xlim([0.1,1e8/1e3])
ylim([0,zUpperLim/1e3])
xline(f_reflex/1e3,"","Reflexion frequency") %Plot vertical line at reflexion frequency
yline(ZRefLinear/1e3,"","Reflexion altitude") %Plot horizontal line at reflexion altitude

fig3=figure(3);
plot(XLinear/1e3,ZLinear/1e3)
%Plot of the coordinates of a linear ray path (both axis in km) FIGURE3
yline(ZRefLinear/1e3,"","Reflexion altitude") %Plot horizontal line at reflexion altitude
xlabel('x [km]')
ylabel('z [km]]')
xlim([0.1,xMaximum/1e3]) %End the X axis at xMaximum value

fig4=figure(4);
%the reflection shall be computed with the difference of the flipped vector Zi
plot(NLX/1e3,NLZ/1e3)
yline(NLZ(breakpoint)/1e3,"","Reflexion altitude") %Plot horizontal line at reflexion altitude
xlabel('x [km]')
ylabel('z [km]]')
xlim([0,NLX(2*breakpoint)/1e3]) %End the X axis at xMaximum value

fig5=figure(5);
%Comparisson Figure
plot(XLinear/1e3,ZLinear/1e3,NLX/1e3,NLZ/1e3)
%Plot of the coordinates of a linear ray path (both axis in km) FIGURE3
yline(ZRefLinear/1e3,"","Reflexion altitude (linear aprox)") %Plot horizontal line at reflexion altitude (linear)
yline(NLZ(breakpoint)/1e3,"","Reflexion altitude") %Plot horizontal line at reflexion altitude
xlabel('x [km]')
ylabel('z [km]]')
xlim([0,(NLX(2*breakpoint)/1e3)+50]) %End the X axis at xMaximum value +50km

%STEP 5
[PF1 Z01 neZ1 zUpperLim1 xMaximum1 f_reflex1 XLinear1 ZLinear1 ZRefLinear1 NLX1 NLZ1 breakpoint1] = IonosphereRayTracing([aa ai ar I0 n0 T+100 f ((alpha*pi)/180)])

fig6=figure(6);
%the reflection shall be computed with the difference of the flipped vector Zi
plot(NLX/1e3,NLZ/1e3,NLX1/1e3,NLZ1/1e3)
yline(NLZ(breakpoint)/1e3,"","Reflexion altitude") %Plot horizontal line at reflexion altitude
yline(NLZ1(breakpoint1)/1e3,"","Reflexion altitude +100K")
xlabel('x [km]')
ylabel('z [km]]')

%Step 7
[PF2 Z02 neZ2 zUpperLim2 xMaximum2 f_reflex2 XLinear2 ZLinear2 ZRefLinear2 NLX2 NLZ2 breakpoint2] = IonosphereRayTracing([aa ai ar I0*10 n0 T+100 f ((alpha*pi)/180)])

fig7=figure(7);
%the reflection shall be computed with the difference of the flipped vector Zi
plot(NLX/1e3,NLZ/1e3,NLX2/1e3,NLZ2/1e3)
yline(NLZ(breakpoint)/1e3,"","Reflexion altitude") %Plot horizontal line at reflexion altitude
yline(NLZ2(breakpoint2)/1e3,"","Reflexion altitude 10 times I0")
xlabel('x [km]')
ylabel('z [km]]')