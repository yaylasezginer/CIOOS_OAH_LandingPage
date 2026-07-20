
function [data_export_nc, stats] = ExportFilesCreate_profiles(savedir, siteID, siteLat, siteLon, sensorDir, loDir, sscDir, params_sensor, params, units, DataURL, loURL, sscURL)

% read data from downloaded files from sensor and model (SalishSeaCast and LiveOcean) data providers 
% Bring datasets to the same temporal resolution and reformat into a single
% matrix. Return matrix as a csv and nc filetype. 
% Note: This fxn is designed for profiling platforms
% Note: update L66 with local directory for nctoolbox
%
% INPUTS:
% loURL: URL for LiveOcean data access/citation page (string)
% sscURL: URL for ssc data access page (string)
% SensorNetworkLog: table with following columns:
%   savedir: folder directory for export files (string)
%   siteID: Name of site (string)
%   siteLat: Lat in degree decimal (scalar)
%   siteLon: Lon in degree decimal (scalar)
%   sensorDir: filepath and filename for sensor data (string)
%   loDir: filepath and filename for liveocean data (string)
%   sscDir: filepath and filename for SalishSeaCast data (string)
%   params_sensor: parameter names given by data provider for variables measured by sensor (cell array)
%   params: standardized names for parameters measured by sensor. May vary from parameter name given
%    by data provider depending on controlled vocabulary usage (cell array). 
%    possible elements of params: time, depth, temp, sal, dO2, DIC, TA, pCO2, pH
%   units: sensor measurement units (cell array equal in size to params)
%    time: possible units: s_1970, s_1950, string format (copy string format
%    exactly), dn (matlab dn)
%    temp: possible units: degC, K
%    sal: possible units: psu, g/kg, ppt
%    dO2: possible units: mmol/m3, umol/kg, umol/L, mL/L, mg/L, %
%    DIC: possible units: mmol/m3, umol/kg, umol/L
%    TA: possible units: mmol/m3, umol/kg, umol/L
%    pCO2: possible units: uatm 
%   DataURL: URL for in-situ data source citation page (string)
%
% OUTPUTS:
% data_export: reformatted data matrix identical to data contents of csv export file. Interpolated values are not included (table)
% data_export_csv: reformatted data matrix identical nc export file
% contents. Includes interpolated values. For CIOOS Pacific internal use
% only. (table)
% metadata: metadata associated with export file data. 
% 
% Data will always be returned in standardized units:
% time: time in UTC
% T: degC
% S: g/kg
% DIC: mmol/m3
% TA:mmol/m3
% dO2: uM
% pH: ()
% pCO2: uatm
% saturation state: ()
% 
% Notes: prior to running, salishseacast data should be compiled into a
% single .nc file


%% ------------------------------------------------------------------------
% Initialize outputs
% -------------------------------------------------------------------------

data_export_nc = struct();

%% ------------------------------------------------------------------------
% DETERMINE FILE TYPES AND READ DATA
% -------------------------------------------------------------------------

% Prepare to read .nc datafiles
addpath('/Users/yaylasezginer/Documents/MATLAB/nctoolbox') 
setup_nctoolbox

% 1. Read sensor data
[~,~,sensorFileType] = fileparts(sensorDir);

switch sensorFileType
    case '.csv'
        sensorData_x = struct(readtable(sensorDir)); 
    case '.mat'
        sensorData_x = load(sensorDir); 
    case '.nc'
        sensorData_nc = ncdataset(sensorDir); 
        % reformat to a struct
        varSensor = sensorData_nc.variables;
        sensorData_x = struct();
        for v = 1:numel(varSensor)
            sensorData_x.(varSensor{v}) = squeeze(sensorData_nc.data(varSensor{v}));
        end
end
% Convert sensor data parameter names from given names to standardized
% vocabulary
for v = 1:numel(params_sensor)
    sensorData.(params{v}) = sensorData_x.(params_sensor{v});
end

% 2. Read LiveOcean data (always in .nc format). 3D data to be reformatted to
% vector from matrix. If site doesn't overlap with LO model domain, return
% empty vectors

varLo = [{'TIC'},{'oxygen'},{'alkalinity'},{'salt'},{'temp'}]; % native LO variable names
stdNames = [{'DIC'},{'dO2'},{'TA'},{'sal'},{'temp'}]; % controlled vocabulary variable names
loData = struct();
if ~isempty(loDir)
    loData_nc = ncdataset(loDir);
    loData.time = loData_nc.data('ocean_time'); % s since 1970
    loData.dn = datenum(1970,1,1,0,0,loData.time);
    loData.depth = -loData_nc.data('z_rho');
    for v = 1:numel(varLo)
        loData.(stdNames{v}) = squeeze(loData_nc.data(varLo{v}));
    end
else 
    loData.dn = [];
    loData.dn = [];
    for v = 1:numel(varLo)
        loData.(stdNames{v}) = [];
    end
end

% 3. Read SSC data (always in .nc format). If site doesn't overlap with LO model domain, return
% empty vectors

varSsc = [{'dissolved_inorganic_carbon'},{'dissolved_oxygen'},{'total_alkalinity'},{'salinity'},{'temperature'}];
sscData = struct();
if ~isempty(sscDir)
    sscData_nc = ncdataset(sscDir);
    sscData.time = sscData_nc.data('time'); % s since 1970
    sscData.dn = datenum(1970,1,1,0,0,sscData.time);
    sscData.depth = sscData_nc.data('depth');
    for v = 1:numel(varSsc)
        sscData.(stdNames{v}) = squeeze(sscData_nc.data(varSsc{v}));
    end
else
    sscData.dn = [];
    sscData.depth = [];
    for v = 1:numel(varSsc)
        sscData.(stdNames{v}) = [];
    end
end

%% ------------------------------------------------------------------------
% Ensure sensor data is in standardized units
% -------------------------------------------------------------------------

addpath /Users/yaylasezginer/Documents/MATLAB/GSW-Matlab-master/Toolbox
addpath /Users/yaylasezginer/Documents/MATLAB/GSW-Matlab-master/Toolbox/library

sensorData.pres = sw_pres(sensorData.depth, siteLat); % units: db
[nSensorR, nSensorC] = size(sensorData.temp);
[nDepthR, nDepthC] = size(sensorData.depth);
if isequal([nDepthR, nDepthC], [nSensorR, nSensorC])
    sensorData.dens = sw_dens(sensorData.sal, sensorData.temp, sensorData.pres); % kg/m3
elseif nDepthR == nSensorR && nDepthC == 1
    P = repmat(sensorData.pres,1,nSensorC);
    sensorData.dens = sw_dens(sensorData.sal, sensorData.temp, P); % kg/m3
elseif nDepthR == nSensorC && nDepthC == 1
    P = repmat(sensorData.pres',nSensorR,1);
    sensorData.dens = sw_dens(sensorData.sal, sensorData.temp, P); % kg/m3   
end

units_out = cell(size(units));

% Time

time_ind = strcmp(params, 'time');
if any(time_ind)
    time_unit = units{time_ind};
    switch time_unit
        case 's_1970'
            sensorData.dn = datenum(1970, 1,1,0,0,sensorData.time);
        case 'dn'
            sensorData.dn = sensorData.time;
        case 's_1950'
            sensorData.dn = datenum(1950, 1,1,0,0,sensorData.time);
        case 'string'
            sensorData.dn = datenum(sensorData.time, 'string'); % replace string with formatIn
    end
    units_out{time_ind} = 'UTC';
end


% Temperature
temp_ind = strcmp(params,'temp');
if any(temp_ind)
    temp_unit = units{temp_ind};
    switch temp_unit
        case 'degC'
            sensorData_std.temp = sensorData.temp;
        case 'K'
            sensorData_std.temp = sensorData.temp - 273.15;
    end
    units_out{temp_ind} = 'degC';
end

% Salinity
sal_ind = strcmp(params,'sal');
if any(sal_ind)
    sal_unit = units{sal_ind};
    switch sal_unit
        case 'psu'
            sensorData_std.sal = gsw_SP_from_SA(sensorData.sal,sensorData.pres,siteLon,siteLat);
        case 'g/kg'
            sensorData_std.sal = sensorData.sal;
        case 'ppt'
            sensorData_std.sal = sensorData.sal;
    end
    units_out{sal_ind} = 'g/kg';
end

% Dissolved O2
dO2_ind = strcmp(params,'dO2');
if any(dO2_ind)
    dO2_unit = units{dO2_ind};
    switch dO2_unit
        case 'mmol/m3'
            sensorData_std.dO2 = sensorData.dO2;
        case 'umol/L'
            sensorData_std.dO2 = sensorData.dO2;
        case 'umol/kg'
            sensorData_std.dO2 = sensorData.dO2 .* 1000 ./ sensorData.dens;
        case 'mL/L'
            sensorData_std.dO2 = 44.6596 .* sensorData.dO2;
        case 'mg/L'
            sensorData_std.dO2 = 31.25 .* sensorData.dO2;
        case 'mg/kg'
            sensorData_std.dO2 = (sensorData.dO2 .* sensorData.dens)./31.998;
        case 'mL/kg'
            sensorData_std.dO2 = (sensorData.dO2 .* sensorData.dens)./22.392;
    end
    units_out{dO2_ind} = 'mmol/m3';
end

% DIC
DIC_ind = strcmp(params, 'DIC');
if any(DIC_ind)
    DIC_unit = units{DIC_ind};
    switch DIC_unit
        case 'mmol/m3'
            sensorData_std.DIC = sensorData.DIC;
        case 'umol/kg'
            sensorData_std.DIC = sensorData.DIC .* 1000 ./ sensorData.dens;
        case 'umol/L'
            sensorData_std.DIC = sensorData.DIC;
    end
    units_out{DIC_ind} = 'mmol/m3';
end

% TA
TA_ind = strcmp(params, 'TA');
if any(TA_ind)
    TA_unit = units{TA_ind};
    switch TA_unit
        case 'mmol/m3'
            sensorData_std.TA = sensorData.TA;
        case 'umol/kg'
            sensorData_std.TA = sensorData.TA .* 1000 ./ sensorData.dens;
        case 'umol/L'
            sensorData_std.TA = sensorData.TA;
    end
    units_out{TA_ind} = 'mmol/m3';
end

disp(['units converted for ' params'])

%% ------------------------------------------------------------------------
% Calculate carbonate system parameters for LiveOcean and SalishSeaCast
% -------------------------------------------------------------------------

% Sal input to CO2SYS in psu
% TA and DIC input to CO2SYS in umol/kgSw
% No matrix input to CO2SYS

% LiveOcean carbonate system 
% Convert Matrix [time x depth] to array

addpath  /Users/yaylasezginer/Documents/MATLAB/CO2-System-Extd-main/main

if ~isempty(loData.dn)
    loData.pres = sw_pres(loData.depth,siteLat); % units: db
    loData.dens = sw_dens(loData.sal, loData.temp, loData.pres); % kg/m3
    loP = reshape(loData.pres,[],1);
    loSigma = reshape(loData.dens,[],1);
    loData.sal_psu = reshape(gsw_SA_from_SP(loData.sal, loData.pres, siteLon, siteLat),[],1);
    [DATA,~,~]=CO2SYS(reshape(loData.TA,[],1)*1000./loSigma,reshape(loData.DIC,[],1).*1000./loSigma, ...
        1,2,loData.sal_psu,reshape(loData.temp,[],1),reshape(loData.temp,[],1),loP,loP,0,0,0,0,1,10,1,2,2);
    % Return data to matrix format (time x depth)
    loData.pCO2 = reshape(DATA(:,4),numel(loData.time),[]);
    loData.pH = reshape(DATA(:,21),numel(loData.time),[]);
    loData.omega_arag = reshape(DATA(:,18),numel(loData.time),[]);
    loData.omega_calc = reshape(DATA(:,17),numel(loData.time),[]);
else
    loData.pCO2 = [];
    loData.pH = [];
    loData.omega_arag = [];
    loData.omega_calc = [];
end

% SSC carbonate system 

if ~isempty(sscData.dn)
    sscData.pres = sw_pres(sscData.depth, siteLat); % units: db
    sscData.pres = repmat(sscData.pres',numel(sscData.time),1);
    sscData.dens= sw_dens(sscData.sal, sscData.temp, sscData.pres); % kg/m3
    sscP = reshape(sscData.pres, [], 1);
    sscSigma = reshape(sscData.dens, [], 1);
    sscData.sal_psu = reshape(gsw_SA_from_SP(sscData.sal, sscData.pres, siteLon, siteLat),[],1);
    [DATA,~,~]=CO2SYS(reshape(sscData.TA,[],1)*1000./sscSigma,reshape(sscData.DIC,[],1).*1000./sscSigma, ...
        1,2,sscData.sal_psu,reshape(sscData.temp,[],1),reshape(sscData.temp,[],1),sscP,sscP,0,0,0,0,1,10,1,2,2);
    % Return data to matrix format (time x depth)
    sscData.pCO2 = reshape(DATA(:,4),[],numel(sscData.depth));
    sscData.pH = reshape(DATA(:,21),[],numel(sscData.depth));
    sscData.omega_arag = reshape(DATA(:,18),[],numel(sscData.depth));
    sscData.omega_calc = reshape(DATA(:,17),[],numel(sscData.depth));
else
    sscData.pCO2 = [];
    sscData.pH = [];
    sscData.omega_arag = [];
    sscData.omega_calc = [];
end

disp('carbonate chemistry calculated for LO and SSC')

%% ------------------------------------------------------------------------
% Preparing data matrices 
% -------------------------------------------------------------------------
 
% Bring model data to matching time and space resolution as sensor data. 
% Use 1D nearest neighbor interpolation to match sensor time stamps (reduce
% model resolution to match sensor time stamps)
% 
% Use indexing to match disperesed model depth bins to nearest depth
% measured by sensor. Depth gaps in model profile filled with NaN.


% Prepare uniform matrix with size = sensor time x depth dimensions

cut_row = all(isnan(sensorData.dn),2); % time gaps in sensor data
cut_col = all(isnan(sensorData.dn),1); % unmeasured depths
fill_dn = min(sensorData.dn(~cut_row,~cut_col),[],2,'omitmissing'); % start time of each cast
binned_depth = floor(min(sensorData.depth(~cut_col))):floor(max(sensorData.depth(~cut_col))); 

nT = numel(fill_dn); % Matrix rows (x) dimension
nZ = numel(binned_depth); % Matrix column (y) dimension 

% Prepare depth (y) indices to match model depths with sensor depth
sscBottom = all(sscData.temp == 0,1)';
sscCut = any([sscBottom, sscData.depth < binned_depth(1)],2); % cut data below local seafloor depth and shallower than sensor range
sscDepth = sscData.depth(~sscCut); 
nZ_SSC = numel(sscDepth);
[depthOffset_SSC, iSSC2sensor] = min(abs(repmat(binned_depth,nZ_SSC,1) - sscDepth),[],2);

loCut = loData.depth(1,:) >= binned_depth(1);
loDepth = loData.depth(1,loCut)'; % convert matrix to Nx1 array
nZ_LO = numel(loDepth);
[depthOffset_LO, iLO2sensor] = min(abs(repmat(binned_depth,nZ_LO,1) - loDepth),[],2);

% Save matrix dimensions data to nc export dataset
data_export_nc.depth_sensor = binned_depth';
data_export_nc.depth_SSC = sscData.depth(~sscBottom);
data_export_nc.depth_LiveOcean = loData.depth;
data_export_nc.time_sensor = (fill_dn - datenum(1970,1,1))./datenum(0,0,0,0,0,1); % convert back to s since 1970
data_export_nc.time_SSC = (sscData.dn - datenum(1970,1,1))./datenum(0,0,0,0,0,1); % convert back to s since 1970
data_export_nc.time_LiveOcean = (loData.dn - datenum(1970,1,1))./datenum(0,0,0,0,0,1); % convert back to s since 1970

for i = 1:numel(params)
    if strcmp(params{i}, 'time') || strcmp(params{i}, 'depth')
        continue
    else

        % initialize 2D arrays in data_export and data_export_nc (time x depth)
        data_export_nc.([params{i} '_SSC_interpolated']) = nan(nT,nZ);
        data_export_nc.([params{i} '_LiveOcean_interpolated']) = nan(nT,nZ);
        
        % Untouched original data
        data_export_nc.([params{i} '_SSC']) = sscData.(params{i})(:,~sscBottom);
        data_export_nc.([params{i} '_LiveOcean']) = loData.(params{i});
        data_export_nc.([params{i} '_sensor']) = sensorData_std.(params{i})(~cut_row,~cut_col);  

        if ~isempty(loData.dn) && ~isempty(loData.(params{i}))
            for z = 1:nZ_LO
                loFill(:,z) = interp1(loData.dn, data_export_nc.([params{i} '_LiveOcean'])(:,z),fill_dn,'nearest',NaN);
            end
            data_export_nc.([params{i} '_LiveOcean_interpolated'])(:,iLO2sensor) = loFill;
        end

        if ~isempty(sscData.dn) && ~isempty(sscData.(params{i}))
            for z = 1:nZ_SSC
                sscFill(:,z) = interp1(sscData.dn, data_export_nc.([params{i} '_SSC'])(:,z),fill_dn,'nearest',NaN);
            end
            data_export_nc.([params{i} '_SSC_interpolated'])(:,iSSC2sensor) = sscFill;
        end


    end

%% ------------------------------------------------------------------------
% Calculating Model Evaluation Stats & Plotting/Saving figures
% -------------------------------------------------------------------------
    obs = data_export_nc.([params{i} '_sensor']);
    stats.LiveOcean.Residuals.(params{i}) = obs - data_export_nc.([params{i} '_LiveOcean_interpolated']);
    stats.SSC.Residuals.(params{i}) = obs - data_export_nc.([params{i} '_SSC_interpolated']);
    nonanLO = isnan(stats.LiveOcean.Residuals.(params{i}));
    nonanSSC = isnan(stats.SSC.Residuals.(params{i}));
    stats.LiveOcean.RMSE.(params{i}) = sqrt(sum((stats.LiveOcean.Residuals.(params{i})(~nonanLO)).^2)./sum(~nonanLO));
    stats.SSC.RMSE.(params{i}) = sqrt(sum((stats.SSC.Residuals.(params{i})(~nonanSSC)).^2)./sum(~nonanSSC));

    sp = figure;
    sp(1) = subplot(2,3,1); 
    imAlpha = ones(size(data_export_nc.([params{i} '_LiveOcean_interpolated'])'));
    imAlpha(isnan(data_export_nc.([params{i} '_LiveOcean_interpolated'])')) = 0;
    imagesc(fill_dn, data_export_nc.depth_sensor, data_export_nc.([params{i} '_LiveOcean_interpolated'])','AlphaData',imAlpha); shading flat; hold on
    ylabel('Depth'); xlabel('time')
    c = colorbar; ylabel(c, [params{i} ' ' units_out{i}])
    datetick('x')
    colormap(cmocean('thermal'))
    title('Live Ocean Interpolated')  
    set(gca,'color',[0.5 0.5 0.5],'FontSize',15,'XLim',[min(fill_dn) max(fill_dn)]);

    sp(2) = subplot(2,3,2);
    imAlpha = ones(size(data_export_nc.([params{i} '_SSC_interpolated'])'));
    imAlpha(isnan(data_export_nc.([params{i} '_SSC_interpolated'])')) = 0;
    imagesc(fill_dn, data_export_nc.depth_sensor, data_export_nc.([params{i} '_SSC_interpolated'])','AlphaData',imAlpha); shading flat; hold on
    ylabel('Depth'); xlabel('time')
    title('SalishSeaCast Interpolated')
    c = colorbar; ylabel(c, [params{i} ' ' units_out{i}])
    datetick('x')
    colormap(cmocean('thermal'))
    set(gca,'color',[0.5 0.5 0.5],'FontSize',15,'XLim',[min(fill_dn) max(fill_dn)])
   
    sp(3) = subplot(2,3,3);
    imagesc(fill_dn, data_export_nc.depth_sensor, data_export_nc.([params{i} '_sensor'])'); shading flat; hold on
    ylabel('Depth'); xlabel('time')
    title('Sensor')
    c = colorbar; ylabel(c, [params{i} ' ' units_out{i}])
    datetick('x')
    colormap(cmocean('thermal'))
    set(gca,'color',[0.5 0.5 0.5],'FontSize',15,'XLim',[min(fill_dn) max(fill_dn)])
  
    sp(4) = subplot(2,3,4);
    imAlpha = ones(size(stats.LiveOcean.Residuals.(params{i})'));
    imAlpha(isnan(stats.LiveOcean.Residuals.(params{i})')) = 0;
    imagesc(fill_dn, data_export_nc.depth_sensor, stats.LiveOcean.Residuals.(params{i})','AlphaData', imAlpha); shading flat; hold on
    ylabel('Depth'); xlabel('time')
    title('Sensor - LiveOcean')
    cmax = max(abs(stats.LiveOcean.Residuals.(params{i})),[],'all','omitnan');
    c = colorbar; ylabel(c, ['\Delta ' units_out{i}]); clim([-cmax cmax])
    datetick('x')
    colormap(cmocean('balance'))
    set(gca,'color',[0.5 0.5 0.5],'FontSize',15,'XLim',[min(fill_dn) max(fill_dn)])

    sp(5) = subplot(2,3,5);
    imAlpha = ones(size(stats.SSC.Residuals.(params{i})'));
    imAlpha(isnan(stats.SSC.Residuals.(params{i})')) = 0;
    imagesc(fill_dn, data_export_nc.depth_sensor,stats.SSC.Residuals.(params{i})','AlphaData',imAlpha); shading flat; hold on
    ylabel('Depth'); xlabel('time')
    title('Sensor - SalishSeaCast')
    cmax = max(abs(stats.SSC.Residuals.(params{i})),[],'all','omitnan');
    c = colorbar; ylabel(c, ['\Delta ' units_out{i}]); clim([-cmax cmax])
    datetick('x')
    colormap(cmocean('balance'))
    set(gca,'color',[0.5 0.5 0.5],'FontSize',15,'XLim',[min(fill_dn) max(fill_dn)])
    
    saveas(gcf,[savedir '/figures/' siteID '_' params{i} '.png'])
end
%% ------------------------------------------------------------------------
% Format metadata & save export files
% -------------------------------------------------------------------------

fileSavePath = [savedir '/' siteID];
[nT_LO, nZ_LO] = size(loData.depth);
nT_SSC = numel(sscData.dn); nZ_SSC = numel(sscData.depth(~sscBottom));

% Formatting for internal .nc files Global attributes
nccreate([fileSavePath '.nc'], 'time_sensor','dimensions',{'t_sensor',nT})
ncwrite([fileSavePath '.nc'], 'time_sensor',data_export_nc.time_sensor)
ncwriteatt([fileSavePath '.nc'],'time_sensor','units','s since 1970,01,01 00:00 (UTC)')

nccreate([fileSavePath '.nc'], 'time_SalishSeaCast','dimensions',{'t_SSC',nT_SSC})
ncwrite([fileSavePath '.nc'], 'time_SalishSeaCast',data_export_nc.time_SSC)
ncwriteatt([fileSavePath '.nc'],'time_SalishSeaCast','units','s since 1970,01,01 00:00 (UTC)')

nccreate([fileSavePath '.nc'], 'time_LiveOcean','dimensions',{'t_LO',nT_LO})
ncwrite([fileSavePath '.nc'], 'time_LiveOcean',data_export_nc.time_LiveOcean)
ncwriteatt([fileSavePath '.nc'],'time_LiveOcean','units','s since 1970,01,01 00:00 (UTC)')

nccreate([fileSavePath '.nc'], 'depth_sensor','dimensions',{'z_sensor',nZ})
ncwrite([fileSavePath '.nc'], 'depth_sensor',data_export_nc.depth_sensor)
ncwriteatt([fileSavePath '.nc'],'depth_sensor','units','m')

nccreate([fileSavePath '.nc'], 'depth_SalishSeaCast','dimensions',{'z_SSC',nZ_SSC})
ncwrite([fileSavePath '.nc'], 'depth_SalishSeaCast',data_export_nc.depth_SSC)
ncwriteatt([fileSavePath '.nc'],'depth_SalishSeaCast','units','m')

nccreate([fileSavePath '.nc'], 'depth_LiveOcean','dimensions',{'t_LO',nT_LO,'z_LO',nZ_LO})
ncwrite([fileSavePath '.nc'], 'depth_LiveOcean',data_export_nc.depth_LiveOcean)
ncwriteatt([fileSavePath '.nc'],'depth_LiveOcean','units','m')

ncwriteatt([fileSavePath '.nc'],'/','Site_ID',siteID);
ncwriteatt([fileSavePath '.nc'],'/','Latitude',siteLat);
ncwriteatt([fileSavePath '.nc'],'/','Longitude',siteLon);
ncwriteatt([fileSavePath '.nc'],'/','Depth','profile');


for i = 1:numel(params)
    if strcmp(params{i}, 'time') || strcmp(params{i}, 'depth')
        continue
    end

    nccreate([fileSavePath '.nc'],[params{i} '_sensor'],'dimensions',{'t_sensor',nT,'z_sensor',nZ},'FillValue','disable')
    ncwrite([fileSavePath '.nc'],[params{i} '_sensor'],data_export_nc.([params{i} '_sensor']))
    ncwriteatt([fileSavePath '.nc'],[params{i} '_sensor'], 'units', units_out{i})
    ncwriteatt([fileSavePath '.nc'],[params{i} '_sensor'], 'source', DataURL)

    nccreate([fileSavePath '.nc'],[params{i} '_SSC'],'dimensions',{'t_SSC',nT_SSC,'z_SSC',nZ_SSC},'FillValue','disable')
    ncwrite([fileSavePath '.nc'],[params{i} '_SSC'],data_export_nc.([params{i} '_SSC']))
    ncwriteatt([fileSavePath '.nc'],[params{i} '_SSC'], 'units', units_out{i})
    ncwriteatt([fileSavePath '.nc'],[params{i} '_SSC'], 'source', sscURL)    

    nccreate([fileSavePath '.nc'],[params{i} '_LiveOcean'],'dimensions',{'t_LO',nT_LO,'z_LO',nZ_LO},'FillValue','disable')
    ncwrite([fileSavePath '.nc'],[params{i} '_LiveOcean'],data_export_nc.([params{i} '_LiveOcean']))
    ncwriteatt([fileSavePath '.nc'],[params{i} '_LiveOcean'], 'units', units_out{i})
    ncwriteatt([fileSavePath '.nc'],[params{i} '_LiveOcean'], 'source', loURL)    

    nccreate([fileSavePath '.nc'],[params{i} '_LiveOcean_resized'],'dimensions',{'t_sensor',nT,'z_sensor',nZ},'FillValue','disable')
    ncwrite([fileSavePath '.nc'],[params{i} '_LiveOcean_resized'],data_export_nc.([params{i} '_LiveOcean_interpolated']))
    ncwriteatt([fileSavePath '.nc'],[params{i} '_LiveOcean_resized'], 'units', units_out{i})
    ncwriteatt([fileSavePath '.nc'],[params{i} '_LiveOcean_resized'], 'source', 'LiveOcean data subsampled to match sensor spatiotemporal resolution')  

    nccreate([fileSavePath '.nc'],[params{i} '_SSC_resized'],'dimensions',{'t_sensor',nT,'z_sensor',nZ},'FillValue','disable')
    ncwrite([fileSavePath '.nc'],[params{i} '_SSC_resized'],data_export_nc.([params{i} '_SSC_interpolated']))
    ncwriteatt([fileSavePath '.nc'],[params{i} '_SSC_resized'], 'units', units_out{i})
    ncwriteatt([fileSavePath '.nc'],[params{i} '_SSC_resized'], 'source', 'SSC data subsampled to match sensor spatiotemporal resolution.')  
end

end