% Modified from copmeasure.m
% Do not fit with the Prieto result, i.e., use all data length for temporal parameters, and use the original sampling frequency
% On June 7, 2008 by Kei
% Original code by Vivian on September 27
% Modified by Kei on September 28
% 'signals' as the input matrix. 1ch-copml, 2ch-copap
% 'trialname' as the name of the trial for output filename
% 'allparameters' as the output vector

% modified by Ray during summer 2017 to only include functionality to
% calculate Prieto's parameters, scroll to the bottom to see the list of
% parameters calculated

% modified by Kei on 170823
% no downsample, input Fs, nfft, fc
% expecting signals has two columns with ML, AP in this order.

% Outputs parameters are in the following order:

function allparameters = copmeasures4(signals, Fs, nfft, fc) 

%frequencies description
% Fsorg = 1000; % sampling frequency of the loaded data
% Fs = 100;   % down sampled frequency
% fc = 5;         %cut off frequency
fspel = 0.1; % frequency for spectrum lower band
fspeu = 5; % frequency for spectrum upper band
% nfft = 2048; %for 100Hz
% nfft = 16384; % for 1000Hz

%*****************************************************
order = 2;          %order of the filter
% [B, A] = butter(order, fc/(Fsorg/2));  %butterworth filter
[B, A] = butter(order, fc/(Fs/2));  %butterworth filter
% T = 70;            %Period of data of the overall file
T = length(signals(:,1))/Fs;            %Period of data of the overall file
z05 = 1.645;        %constant for the z statistic at the 95% confidence level
F05 = 3;            %F statistic at 95% confidence level for a bivariate distribution for large sample size

%assigning data to variables
%cutting first 5 seconds and last 5 seconds from 70 s trial
% xCOP  = signals(5001:65000,1); %converting units
% yCOP = signals(5001:65000,2); 
xCOP  = signals(:,1);
yCOP = signals(:,2); 

%filtering signals
%xCOPf = filtfilt(B, A, xCOP);
%yCOPf = filtfilt(B, A, yCOP);
%no filtering because COP data was already filtered by Angela -Ray
xCOPf = xCOP;
yCOPf = yCOP;

%downsample to match Prieto's analysis
% xCOPfd = downsample(xCOPf, Fsorg/Fs);
% yCOPfd = downsample(yCOPf, Fsorg/Fs);
xCOPfd = xCOPf;
yCOPfd = yCOPf;

%using unfiltered, zero-mean data
% APnf = yCOP;
% MLnf = xCOP;
% RDnf = sqrt(APnf.^2+MLnf.^2);
%using unfiltered, zero-mean data
AP = yCOPfd - mean(yCOPfd);
ML = xCOPfd - mean(xCOPfd);
RD = sqrt(AP.^2+ML.^2);

meanCOPy=mean(yCOPfd);
meanCOPx=mean(xCOPfd);
%************************************************************
%Temporal Domain Measures
%Mean distance  
MDISTrd = mean(RD);
MDISTap = mean(abs(AP));
MDISTml = mean(abs(ML));

%RMS distance
sdrd = sqrt(mean(RD.^2));
sdap = sqrt(mean(AP.^2));     %also equal RDISTap
sdml = sqrt(mean(ML.^2));     %also equal RDISTml

%Total Excursion
TOTEXrd = sum(sqrt(diff(AP).^2 + diff(ML).^2));
TOTEXap = sum(abs(diff(AP)));
TOTEXml = sum(abs(diff(ML)));

%Mean velocity
MVELOrd = TOTEXrd/T;
MVELOap = TOTEXap/T;
MVELOml = TOTEXml/T;

%Range
RANGErd = max(RD)-min(RD);
RANGEap = max(AP)-min(AP);
RANGEml = max(ML)-min(ML);

%95% confidence circle AREA-CC
%assuming that distances are normally distributed; CC = %confidence circle
srd = sqrt(sdrd.^2 - MDISTrd.^2);
AREACC = pi * (MDISTrd + z05 * srd)^2;

%95% confidence ellipse - AREA-CE
sdapml = mean(AP.*ML);
AREACE = 2 * pi * F05 * sqrt(sdap^2 * sdml^2 - sdapml^2);

%Sway area AREA-SW
AREASWt = 0;
for i = 1:length(AP)-1
    AREASWt = AREASWt + (abs(AP(i+1)*ML(i) - AP(i)*ML(i+1) ));
end
AREASW = AREASWt/(2*T);

%Mean frequency
MFREQrd = MVELOrd / (2 * pi * MDISTrd);
MFREQap = MVELOap / (4 * sqrt(2) * MDISTap);
MFREQml = MVELOml / (4 * sqrt(2) * MDISTml);

%************************************************************
%Frequency Domain Measures

N = length(AP);
%nfft = max(256, 2^nextpow2(N));
%    nfft = 2048;
noverlap = 0;
%seg = N / nfft;
% [Pap, f] = pwelch(AP, 2048, noverlap, nfft, Fs);
% [Pml, f] = pwelch(ML, 2048, noverlap, nfft, Fs);
% [Prd, f] = pwelch(RD, 2048, noverlap, nfft, Fs);
% [Pap, f] = pwelch(AP, nfft, noverlap, nfft, Fs);
% [Pml, f] = pwelch(ML, nfft, noverlap, nfft, Fs);
% [Prd, f] = pwelch(RD, nfft, noverlap, nfft, Fs);
% [Pap, f] = pwelch(AP, 2000, noverlap, nfft, Fs);
% [Pml, f] = pwelch(ML, 2000, noverlap, nfft, Fs);
% [Prd, f] = pwelch(RD, 2000, noverlap, nfft, Fs);
[Pap, f] = pwelch(AP, nfft, noverlap, nfft, Fs);
[Pml, f] = pwelch(ML, nfft, noverlap, nfft, Fs);
[Prd, f] = pwelch(RD, nfft, noverlap, nfft, Fs);

%find the index where fspel < f < fspeu
temp1 = find(f >= fspel);
temp2 = find(f >= fspeu);
ideltaf = linspace(temp1(1), temp2(1), temp2(1)-1);
% ideltaf = linspace(2, temp2(1), temp2(1)-1); % eliminate only zero Hz
% if Fs == 100
%     temp = find(f >= fc);
%     % In case Fs=100, line 4 = 0.15Hz
%     ideltaf = linspace(4, temp(1), temp(1)-1);
% elseif Fs == 1000
%     temp = find(f >= fc);
%     % In case Fs=1000, line 3 = 0.12Hz
%     ideltaf = linspace(3, temp(1), temp(1)-1);
% end
    
%define the frequency vector for calculating spectral moments
freq = f(ideltaf(1): ideltaf(end));

%define the power spectral density vectors up to fc
Pap1 = Pap(ideltaf(1) : ideltaf(end));
Pml1 = Pml(ideltaf(1) : ideltaf(end));
Prd1 = Prd(ideltaf(1) : ideltaf(end));

%zero moment
u0ap = sum( Pap1 );
u0ml = sum( Pml1 );
u0rd = sum( Prd1 );

%first moment
u1ap = sum( freq.*Pap1 );
u1ml = sum( freq.*Pml1 );
u1rd = sum( freq.*Prd1 );

%second moment
u2ap = sum( (freq.^2).*Pap1 );
u2ml = sum( (freq.^2).*Pml1 );
u2rd = sum( (freq.^2).*Prd1 );

%     %zero moment
%     u0ap(k) = sum( Pap1 ) * f(2);
%     u0ml(k) = sum( Pml1 ) * f(2);
%     u0rd(k) = sum( Prd1 ) * f(2);
%  
%     %first moment
%     u1ap = sum( freq.*Pap1 ) * f(2);
%     u1ml = sum( freq.*Pml1 ) * f(2);
%     u1rd = sum( freq.*Prd1 ) * f(2);
%  
%     %second moment
%     u2ap = sum( (freq.^2).*Pap1 ) * f(2);
%     u2ml = sum( (freq.^2).*Pml1 ) * f(2);
%     u2rd = sum( (freq.^2).*Prd1 ) * f(2);
% 
%centroidal frequency
CFREQap = sqrt(u2ap/u0ap); 
CFREQml = sqrt(u2ml/u0ml); 
CFREQrd = sqrt(u2rd/u0rd); 

%frequency dispersion
FREQDap = sqrt( 1 - (u1ap^2 / (u0ap * u2ap)));
FREQDml = sqrt( 1 - (u1ml^2 / (u0ml * u2ml)));
FREQDrd = sqrt( 1 - (u1rd^2 / (u0rd * u2rd)));

%determining the 50% and 95% power
P50ap = 0.5*u0ap; 
P50ml = 0.5*u0ml;
P50rd = 0.5*u0rd;
P95ap = 0.95*u0ap;
P95ml = 0.95*u0ml;
P95rd = 0.95*u0rd;    

%temporary markers for storing the index at which the 50% and 95% power
%frequencies occur
marker50ap = 0;
marker50ml = 0;
marker50rd = 0;
marker95ap = 0;
marker95ml = 0;
marker95rd = 0;

for i = 1: length(freq)
    if marker50ap == 0
        if sum(Pap1(1:i)) >= P50ap
%             if sum(Pap1(1:i))*freq(1) >= P50ap
            marker50ap = i;
        end
    end
    if marker95ap == 0
        if sum(Pap1(1:i)) >= P95ap
%             if sum(Pap1(1:i))*freq(1) >= P95ap
            marker95ap = i;
        end
    end
    if marker50ml == 0
        if sum(Pml1(1:i)) >= P50ml
%             if sum(Pml1(1:i))*freq(1) >= P50ml
            marker50ml = i;
        end
    end
    if marker95ml == 0
        if sum(Pml1(1:i))>= P95ml
%             if sum(Pml1(1:i))*freq(1) >= P95ml
            marker95ml = i;
        end
    end
    if marker50rd == 0
        if sum(Prd1(1:i)) >= P50rd
%             if sum(Prd1(1:i))*freq(1) >= P50rd
            marker50rd = i;
        end
    end
    if marker95rd == 0
        if sum(Prd1(1:i)) >= P95rd
%             if sum(Prd1(1:i))*freq(1) >= P95rd
            marker95rd = i;
        end
    end        
end

%50% power frequency
PF50ap = freq(marker50ap);
PF50ml = freq(marker50ml);
PF50rd = freq(marker50rd);

%95% power frequency
PF95ap = freq(marker95ap);
PF95ml = freq(marker95ml);
PF95rd = freq(marker95rd);


%**********************************************************************
%Stabilogram Diffusion Analysis

N = length(AP);

%calculating the stabilogram
for m = 1 : Fs*10+1
% for m = 1 : 1001
    DISPrd = 0;
    DISPap = 0;
    DISPml = 0;
    count = 0;
    for i = 1: (N-(m-1))
        DISPrd =  (RD(i+(m-1))-RD(i))^2 + DISPrd ;
        DISPap =  (AP(i+(m-1))-AP(i))^2 + DISPap ;
        DISPml =  (ML(i+(m-1))-ML(i))^2 + DISPml ;
        count = count + 1;
    end
    RDISTrd(m) = 1/count * DISPrd;
    RDISTap(m) = 1/count * DISPap;
    RDISTml(m) = 1/count * DISPml;
end

t = linspace(0, 10, Fs*10+1);
% t = linspace(0, 10, 1001);

%***************
%defining the regions for linear fit
%short term region
a = find(t == 0.5);
pstrd = polyfit(t(1:a), RDISTrd(1:a), 1);
pstap = polyfit(t(1:a), RDISTap(1:a), 1);
pstml = polyfit(t(1:a), RDISTml(1:a), 1);

%long term region
b = find(t == 2);
pltrd = polyfit(t(b:end), RDISTrd(b:end), 1);
pltap = polyfit(t(b:end), RDISTap(b:end), 1);
pltml = polyfit(t(b:end), RDISTml(b:end), 1);

%calculating the coefficients for the linear fit
%short term fit
fstrd = polyval(pstrd, t);
fstap = polyval(pstap, t);
fstml = polyval(pstml, t);

%long term fit
fltrd = polyval(pltrd, t);
fltap = polyval(pltap, t);
fltml = polyval(pltml, t);

%calculating the D values
Dsrd = 0.5*pstrd(1);
Dsap = 0.5*pstap(1);
Dsml = 0.5*pstml(1);

Dlrd = 0.5*pltrd(1);
Dlap = 0.5*pltap(1);
Dlml = 0.5*pltml(1);

%calculating critical points
Dcrd = -(pstrd(2)-pltrd(2))/(pstrd(1)-pltrd(1));    %critical point in terms of time
Dcap = -(pstap(2)-pltap(2))/(pstap(1)-pltap(1));    %critical point in terms of time
Dcml = -(pstml(2)-pltml(2))/(pstml(1)-pltml(1));    %critical point in terms of time
Dintrd = polyval(pstrd, Dcrd);    %critical point in terms of "y"
Dintap = polyval(pstap, Dcap);    %critical point in terms of "y"
Dintml = polyval(pstml, Dcml);    %critical point in terms of "y"

%***************
%calculate log-log plots
tlog = log10(t(2:end));
RDISTrdlog = log10(RDISTrd(2:end));
RDISTaplog = log10(RDISTap(2:end));
RDISTmllog = log10(RDISTml(2:end));

%short term region
c1 = find(tlog == log10(0.01));
c = find(tlog == log10(0.5));
pstrdlog = polyfit(tlog(c1:c), RDISTrdlog(c1:c), 1);
pstaplog = polyfit(tlog(c1:c), RDISTaplog(c1:c), 1);
pstmllog = polyfit(tlog(c1:c), RDISTmllog(c1:c), 1);

%long term region
d = find(tlog == log10(2));    
pltrdlog = polyfit(tlog(d:end), RDISTrdlog(d:end), 1);
pltaplog = polyfit(tlog(d:end), RDISTaplog(d:end), 1);
pltmllog = polyfit(tlog(d:end), RDISTmllog(d:end), 1);

%short term fit


fstrdlog = polyval(pstrdlog, tlog);
fstaplog = polyval(pstaplog, tlog);
fstmllog = polyval(pstmllog, tlog);

%long term fit
fltrdlog = polyval(pltrdlog, tlog);
fltaplog = polyval(pltaplog, tlog);
fltmllog = polyval(pltmllog, tlog);


%calculating H values
Hsrd = 0.5*pstrdlog(1);
Hsap = 0.5*pstaplog(1);
Hsml = 0.5*pstmllog(1);

Hlrd = 0.5*pltrdlog(1);
Hlap = 0.5*pltaplog(1);
Hlml = 0.5*pltmllog(1);

%calculating critical points
Hcrd = -(pstrdlog(2)-pltrdlog(2))/(pstrdlog(1)-pltrdlog(1));    %critical point in terms of time
Hcap = -(pstaplog(2)-pltaplog(2))/(pstaplog(1)-pltaplog(1));    %critical point in terms of time
Hcml = -(pstmllog(2)-pltmllog(2))/(pstmllog(1)-pltmllog(1));    %critical point in terms of time


Hintrd = polyval(pstrdlog, Hcrd);    %critical point in terms of "y"
Hintap = polyval(pstaplog, Hcap);    %critical point in terms of "y"
Hintml = polyval(pstmllog, Hcml);    %critical point in terms of "y"


%***************
% %plotting stabilograms
% figure
% subplot(3,1,1)
% plot(t, RDISTrd, 'b', t, fstrd, 'k--', t, fltrd, 'k--')
% hold on
% plot(Dcrd, Dintrd, 'r*')
% title('RDISTrd')
% subplot(3,1,2)
% plot(t, RDISTap, 'b', t, fstap, 'k--', t, fltap, 'k--')
% hold on
% plot(Dcap, Dintap, 'r*')
% title('RDISTap')
% subplot(3,1,3)
% plot(t, RDISTml, 'b', t, fstml, 'k--', t, fltml, 'k--')
% hold on
% plot(Dcml, Dintml, 'r*')
% title('RDISTml')    
% 
% figure
% subplot(3,1,1)
% plot(tlog, RDISTrdlog, 'b', tlog, fstrdlog, 'k--', tlog, fltrdlog, 'k--')
% hold on
% plot(Hcrd, Hintrd, 'r*')
% title('RDISTrd log')
% subplot(3,1,2)
% plot(tlog, RDISTaplog, 'b', tlog, fstaplog, 'k--', tlog, fltaplog, 'k--')
% hold on
% plot(Hcap, Hintap, 'r*')
% title('RDISTap log')    
% subplot(3,1,3)
% plot(tlog, RDISTmllog, 'b', tlog, fstmllog, 'k--', tlog, fltmllog, 'k--')
% hold on
% plot(Hcml, Hintml, 'r*')
% title('RDISTml log')    
% 
% pause
% close all

%***************
%Output stabilogram difusion function for each file

%{
Commented out because don't want outfile at this moment -Ray
_____________________________________________________________
ofilename2=[trial,'sdfd','.rlt']; % Outfile name
OutputMeasure2=[t' RDISTrd' fstrd' fltrd' RDISTap' fstap' fltap' RDISTml' fstml' fltml'];
save(ofilename2,'OutputMeasure2','-ASCII');
ofilename3=[trial,'sdfh','.rlt']; % Outfile name
OutputMeasure3=[tlog' RDISTrdlog' fstrdlog' fltrdlog' RDISTaplog' fstaplog' fltaplog' RDISTmllog' fstmllog' fltmllog'];
save(ofilename3,'OutputMeasure3','-ASCII');
%}

%**********************************************************************
%Output all parameters
% MDISTrd MDISTap MDISTml sdrd sdap sdml MVELOrd MVELOap MVELOml TOTEXrd TOTEXap TOTEXml RANGErd RANGEap RANGEml AREACC AREACE AREASW MFREQrd MFREQap MFREQml u0rd u0ap u0ml PF50rd PF50ap PF50ml PF95rd PF95ap PF95ml CFREQrd CFREQap CFREQml FREQDrd FREQDap FREQDml Dcrd Dcap Dcml Dsrd Dsap Dsml Dlrd Dlap Dlml Hcrd Hcap Hcml Hsrd Hsap Hsml Hlrd Hlap Hlml meanCOPy meanCOPx
allparameters=[ MDISTrd MDISTap MDISTml sdrd sdap sdml ...
                MVELOrd MVELOap MVELOml TOTEXrd TOTEXap TOTEXml ...
                RANGErd RANGEap RANGEml AREACC AREACE AREASW ...
                MFREQrd MFREQap MFREQml u0rd u0ap u0ml ...
                PF50rd PF50ap PF50ml PF95rd PF95ap PF95ml ...
                CFREQrd CFREQap CFREQml FREQDrd FREQDap FREQDml...
                Dcrd Dcap Dcml Dsrd Dsap Dsml Dlrd Dlap Dlml Hcrd Hcap Hcml Hsrd Hsap Hsml Hlrd Hlap Hlml...
                meanCOPy meanCOPx];
