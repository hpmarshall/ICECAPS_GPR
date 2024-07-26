%folder='/Users/hpmarshall/HP_DRIVE/ICECAPS_GPR/2024/level0/5Ghz/';
folder='/Users/hpmarshall/HP_DRIVE/ICECAPS_GPR/2024/rawgprdatafromlatejunemelt/melt_event_late_june_5Ghz/'
list=dir([folder,'*.log']); list={list.name}'; %list=list(end-1:end);

%% read raw data files
for i=1:length(list)
    filename=[folder,list{i}];
    try
        [datar, config, switch_config] = cilantro_ancho_logger_parse_RLH(filename);
    catch
        disp(['error on filename:' filename])
        f(i)=NaN;
        DM(1:size(d2,2),i)=NaN*zeros(size(d2'));
    end
        %sample delay = 3.867e-09 !! offset of 3.8 ns
    f(i)=config.SamplesPerSecond;
    d2=f_DAC2Volt(mean(datar,1),config); % convert DAC to volts
    DM(1:size(d2,2),i)=d2';
    i
end
%% dates
for i=1:length(list)
    data.dates(i,:)=list{i}(21:32);
end
data.dates=datetime(str2num(data.dates(:,1:2))+2000,str2num(data.dates(:,3:4)),...
    str2num(data.dates(:,5:6)),str2num(data.dates(:,8:9)),str2num(data.dates(:,11:12)),zeros(length(list),1));
data.list=list;
data.f=f;

%% lets plot the data at this stage
figure(1);clf
imagesc(DM,[0.5 0.6]); colorbar

%% sampling frequency correction
tot_time=data.f.^(-1).*size(DM,1);
e1=e_snowdry(350,5e9,-10); % estimate of dielectric constant
v=3e8/sqrt(real(e1)); % velocity in snow
dsnow=v*tot_time/2;
d0=config.SampleDelay*v/2;  % sample delay offset
DM2=NaN(size(DM));
for t=1:size(DM,2)
    dscale(:,t)=linspace(d0,dsnow(t),size(DM,1));
end
max_d=max(dscale(end,:));
% interpolate to correct for different sampling frequencies
dnew=linspace(d0,max_d,size(DM,1)); %%%%9.9banner 11.5bogus 9.65bogus2
for i=1:size(DM,2) %loop over all traces
    DM2(:,i)=interp1(dscale(:,i),DM(:,i),dnew);
    i
end
%%
figure(2);clf
imagesc(data.dates,dnew,DM2,[0.5 0.6]); colorbar


%%
Ix=find(dnew>3.74); % Use 4.85 if using velocity in air);
DM2=DM2(1:Ix(1),:);
dnew=dnew(1:Ix(1))

%% bandpass filter
dt=median(diff(dnew)*2./v); % get sampling interval
wavelength=7e9;
fmin=wavelength-3*10^9;
fmax=wavelength+3*10^9;
DM3=filter_bandpass(DM2,dt,fmin,fmax);
%% remove first 3 traces (bad data)
data.dates=data.dates(4:end);
DM3=DM3(:,4:end);
DM3env=envelope(DM3,10,'rms');
%%
figure(3);clf
subplot(1,2,1)
imagesc(data.dates,dnew,DM3,[0 0.1]); colorbar
set(gca,'LineWidth',2,'FontSize',14,'FontWeight','bold')
xlabel('date')
ylabel('distance in snow [m]')
subplot(1,2,2)
imagesc(data.dates,dnew,DM3env,[0 0.1]); colorbar
set(gca,'LineWidth',2,'FontSize',14,'FontWeight','bold')
xlabel('date')
ylabel('distance in snow [m]')



% %% try median subtraction before envelope
% bg=nanmedian(DM3,2)*ones(1,size(DM3,2)); % median subtraction
% DM4=DM3-bg2; % background subtraction
% DM4(DM4<0)=0;


%%
% 
% %%
% figure(4);clf
% imagesc(data.dates,dnew,DM4); colorbar

% %% background removal
% bg=nanmedian(DM4,2)*ones(1,size(DM4,2)); % median subtraction
% bg2=prctile(DM4,5,2)*ones(1,size(DM4,2)); % lower 5% subtraction
% DM5=DM4-bg2; % background subtraction
% DM5(DM5<0)=0;
% %% distance correction/range gain
% DM6=DM5.^2.*(dnew'*ones(1,size(DM5,2))).^2;
% 
% %%
% figure(4);clf
% imagesc(data.dates,dnew,10*log10(DM5)); colorbar
% 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% *dist^2 (no calibration applied)
% d=data.d_vert/cosd(40);
% for ch=1:4
%     data.dcorr(:,:,ch)=data.sfcorr2(:,:,ch).^2.*((d'*ones(1,size(data.sfcorr,2))).^2);
% end
% 


% 
% %% plot figure
% % ground surface = 28ns --> 28e9*3e8*cosd(40)/2=3.2m
% % insitu.SD=filloutliers(insitu.SD,'linear','movmedian',20);
% % t_snotel=(3.2-insitu.SD/100)*2/cosd(40)/3e8/1e-9;
% clim=[-70 -30; -70 -30; -75 -35; -75 -35];
% %clim=[-25 -10; -25 -10; -30 -15; -30 -15];
% cmap=colormap((cbrewer('seq','GnBu',64,'spline'))); cmap=cmap.^1.5;%cmap(1:10,:)=[]; %cmap1=fliplr(cmap1);
% figure('Position',[0 0 1000 600]);
% tickss=datetime(datenum(data.dates(1)):10:datenum(data.dates(end)),'ConvertFrom','datenum');%tickss=dateshift(tickss, 'dayofweek','Monday','previous');
% for ch=1:4
%     subplot(2,2,ch)
%     imagesc(datenum(data.dates),data.t_new*1e9,10*log10(data.dcorr(:,:,ch).*data.K(ch))); %ipv data.d*cosd(40).*data.K(ch))
%     colorbar; colormap(cmap);
%     caxis(clim(ch,:));
%     title(data.channels{ch})
%     xticks(datenum(tickss))
%     datetick('x','dd-mmm','keeplimits','keepticks')
%     xtickangle(40)
%     set(gca,'TickDir','out');
%     ylabel('t (ns)')
%     %hold on; plot(datenum(insitu.Date),t_snotel,'k','LineWidth',1)
%     %plot([data.dates(1) data.dates(end)],[3,3],'w','LineWidth',2);%plot([1,121],[670,670]);
% end
% 

%print('BogusCal_2122.png','-dpng','-r500')

% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clim=[-80 -30; -80 -30; -35 -15; -35 -15];
% figure('Position',[0 0 1000 600]);
% for ch=1:2
%     subplot(2,1,ch) %VV
%     %if ch==2, ch=3; end
%     imagesc(datenum(data.dates),data.d_vert,10*log10(data.dcorr(:,:,ch)));%.*data.K(ch)));
%     colorbar; colormap(cmap); caxis(clim(ch,:));
%     title(data.channels{ch})
%     xticks(datenum(tickss)); xtickangle(40)
%     datetick('x','dd-mmm','keeplimits','keepticks')
%     set(gca,'TickDir','out');
%     ylabel('t (ns)')
%     xlim([datenum(data.dates(1)) datenum(data.dates(end))])
%     ylim([1.5 8])
% end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Antenna pattern measurement
% % Location of sphere: sample st-nd
% st=500;
% nd=580;
% az=PanTiltData.Azimuth;
% for ch = 1:4
% data.int(ch,:)=trapz( data.dcorr(st:nd,:,ch) ) ;
% end
%
% meanInt=mean(mean(data.int)); %to calc dBi
% meanInt=max(max(data.int));
%
% figure; hold on;
% for ch=1:4
% plot(10*log10(data.int(ch,:)/meanInt));
% end
% legend(data.channels)
% xlabel('angle'); ylabel('dBi'); xlim([-180 180])
% figure; polarpattern(az(362:end),10*log10(data.int(1,(362:end))/meanInt),az(362:end),10*log10(data.int(2,(362:end))/meanInt))
% %print('tower_asc_desc.png','-dpng')

