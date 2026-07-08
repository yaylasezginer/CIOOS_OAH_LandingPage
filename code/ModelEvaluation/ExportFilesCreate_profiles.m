
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
% Bring data to consistent depth and temporal resolution avoiding gaps
% -------------------------------------------------------------------------

% For 2D data (time x depth), bin/aggregate sensor, LiveOcean, and SSC
% onto the common compiled_dn x common_depth grid. Assume variables:
% - sensorData_std.depth (vector), loData.depth, sscData.depth
% - sensorData_std.<param> is [nt_sensor x nd_sensor], loData.<param> is [nt_lo x nd_lo], sscData.<param> is [nt_ssc x nd_ssc]
% Build common depth vector as union of depths

max_depth = zeros(3,1);
if ~isempty(sscData.depth)
    max_depth(1) = max(sscData.depth(sscData.temp(1,:) ~=0 ));
end
if ~isempty(loData.depth)
    max_depth(2) = max(loData.depth,[],'all');
end
max_depth(3) = max(sensorData.depth,[],'all');
binned_depth = 0:floor(max(max_depth)); 

cut_row = all(isnan(sensorData.dn),2); % time gaps in sensor data
cut_col = all(isnan(sensorData.dn),1); % unmeasured depths
sensorData.castStart_dn = min(sensorData.dn(~cut_row,~cut_col),[],2,'omitmissing');
compiled_dn = unique([sscData.dn; loData.dn; sensorData.castStart_dn]);
nT = numel(compiled_dn);
nZ = numel(binned_depth);

data_export_nc.depth = binned_depth;
data_export_nc.time = (compiled_dn - datenum(1970,1,1))./datenum(0,0,0,0,0,1); % convert back to s since 1970

% Adjustable interpolation thresholds 

maxGapDays = 36/24; %  1 hr = 1/24 day (adjustable)
maxMissing = 0.3; % (1-maxMissing) = maximum allowable fraction of data to be interpolated

for i = 1:numel(params)
    if strcmp(params{i}, 'time') || strcmp(params{i}, 'depth')
        continue
    else
        % initialize 2D arrays in data_export and data_export_nc (time x depth)
        data_export_nc.([params{i} '_sensor']) = nan(nT,nZ);
        data_export_nc.([params{i} '_SSC']) = nan(nT,nZ);
        data_export_nc.([params{i} '_LiveOcean']) = nan(nT,nZ);

        % 
        sensorData_std_trimmed.(params{i}) = sensorData_std.(params{i})(~cut_row,~cut_col);

        % Perform 2D binning by nearest-time matching then depth averaging

        % For high-resolution sensor data: map sensor times to compiled times, then for each compiled time
        % aggregate sensor depths into compiled_depth by averaging data within depth bins.
        if ~isempty(sensorData.dn) && ~isempty(sensorData_std.(params{i}))
            [~,iSensor, iCompiledxSensor] = intersect(sensorData.castStart_dn, compiled_dn);
            % For each matching time, bin sensor depth profiles to compiled_depth
            for k = 1:numel(iCompiledxSensor)
                tidx = iCompiledxSensor(k);
                srcRow = iSensor(k);
                data_export_nc.([params{i} '_sensor'])(tidx,:) = interp1(sensorData.depth(~cut_col), sensorData_std_trimmed.(params{i})(srcRow,:),binned_depth,'linear',NaN);
            end
        end

        % Fill LiveOcean and SSC by matching their time indices to compiled_dn and interpolating/extrapolating in depth
        if ~isempty(loData.dn) && ~isempty(loData.(params{i}))
            [~,iLO, iCompiledxLO] = intersect(loData.dn, compiled_dn);
            for k = 1:numel(iCompiledxLO)
                tidx = iCompiledxLO(k);
                srcRow = iLO(k);
                data_export_nc.([params{i} '_LiveOcean'])(tidx,:) = interp1(loData.depth(srcRow,:), loData.(params{i})(srcRow,:), binned_depth, 'linear', NaN);
            end
            data_export_nc.([params{i} '_LiveOcean_interpolated']) = data_export_nc.([params{i} '_LiveOcean']);
            for z = 1:nZ
                fill_ind = ~isnan(data_export_nc.([params{i} '_sensor'])(:,z));
                fill_dn = compiled_dn(fill_ind);
                data_export_nc.([params{i} '_LiveOcean_interpolated'])(fill_ind,z) = interp1(compiled_dn(iCompiledxLO),data_export_nc.([params{i} '_LiveOcean'])(iCompiledxLO,z),fill_dn,'linear');
            end
        end

        if ~isempty(sscData.dn) && ~isempty(sscData.(params{i}))
            [~,iSSC, iCompiledxSSC] = intersect(sscData.dn, compiled_dn);
            sscData.(params{i})(sscData.(params{i}) == 0) = nan;
            for k = 1:numel(iCompiledxSSC)
                tidx = iCompiledxSSC(k);
                srcRow = iSSC(k);
                data_export_nc.([params{i} '_SSC'])(tidx,:) = interp1(sscData.depth, sscData.(params{i})(srcRow,:), binned_depth, 'nearest', NaN);
            end
            data_export_nc.([params{i} '_SSC_interpolated']) = data_export_nc.([params{i} '_SSC']);
            for z = 1:nZ
                fill_ind = ~isnan(data_export_nc.([params{i} '_sensor'])(:,z));
                fill_dn = compiled_dn(fill_ind);
                data_export_nc.([params{i} '_SSC_interpolated'])(fill_ind,z) = interp1(compiled_dn(iCompiledxSSC),data_export_nc.([params{i} '_SSC'])(iCompiledxSSC,z),fill_dn,'linear');
            end
        end

        % skip_col = all(isnan(data_export_nc.temp_sensor),1);
        % skipped = 0; % Keep track of skipped iterations to make sure the number matches the expected skip_col
        % for z = 1:nZ
        %     valid = find(~isnan(col));
        %     if skip_col(z) & numel(valid) > 
        %         skipped = skipped + 1; 
        %         continue
        %     end
        %     col = data_export_nc.([params{i} '_sensor'])(:,z);
        % 
        %     gapInd = find(diff(compiled_dn(valid)) >= maxGapDays);
        %     gapInd = [gapInd; numel(valid)];
        %     ind = valid(1):valid(gapInd(1)-1);
        %     for chunk = 1:numel(gapInd)-1
        %         if numel(ind) > 3 && sum(~isnan(col(ind))) > maxMissing*numel(ind)
        %             data_export_nc.([params{i} '_sensor'])(ind,z) = fillmissing(col(ind),'linear');
        %         end
        %         if chunk < numel(gapInd)
        %             ind = valid(gapInd(chunk)+1):valid(gapInd(chunk+1)-1);
        %         end
        %     end
        % end
        % disp(['skipped ' num2str(skipped) ' sensor interpolation iterations'])
    end



%% ------------------------------------------------------------------------
% Calculating Model Evaluation Stats & Plotting/Saving figures
% -------------------------------------------------------------------------
    obs = data_export_nc.([params{i} '_sensor']);
    obsNaN = all(isnan(data_export_nc.([params{i} '_sensor'])),2);
    stats.LiveOcean.Residuals.(params{i}) = obs(~obsNaN,:) - data_export_nc.([params{i} '_LiveOcean_interpolated'])(~obsNaN,:);
    stats.SSC.Residuals.(params{i}) = obs(~obsNaN,:) - data_export_nc.([params{i} '_SSC_interpolated'])(~obsNaN,:);
    nonanLO = isnan(stats.LiveOcean.Residuals.(params{i}));
    nonanSSC = isnan(stats.SSC.Residuals.(params{i}));
    stats.LiveOcean.RMSE.(params{i}) = sqrt(sum((stats.LiveOcean.Residuals.(params{i})(~nonanLO)).^2)./sum(~nonanLO));
    stats.SSC.RMSE.(params{i}) = sqrt(sum((stats.SSC.Residuals.(params{i})(~nonanSSC)).^2)./sum(~nonanSSC));

    sp = figure;
    sp(1) = subplot(2,3,1);
    loNaN = all(isnan(data_export_nc.temp_LiveOcean),2);
    pcolor(compiled_dn(~loNaN), -binned_depth, data_export_nc.([params{i} '_LiveOcean'])(~loNaN,:)'); shading flat; hold on
    ylabel('Depth'); xlabel('time')
    c = colorbar; ylabel(c, [params{i} ' ' units_out{i}])
    datetick('x')
    colormap(cmocean('thermal'))
    title('Live Ocean')
    set(gca,'FontSize',15,'XLim',[min(compiled_dn) max(compiled_dn)])

    sp(2) = subplot(2,3,2);
    pcolor(compiled_dn, -binned_depth, data_export_nc.([params{i} '_SSC'])'); shading flat; hold on
    ylabel('Depth'); xlabel('time')
    title('SalishSeaCast')
    c = colorbar; ylabel(c, [params{i} ' ' units_out{i}])
    datetick('x')
    colormap(cmocean('thermal'))
    set(gca,'FontSize',15,'XLim',[min(compiled_dn) max(compiled_dn)])
   
    sp(3) = subplot(2,3,3);
    pcolor(compiled_dn(~obsNaN), -binned_depth, data_export_nc.([params{i} '_sensor'])(~obsNaN,:)'); shading flat; hold on
    ylabel('Depth'); xlabel('time')
    title('Sensor')
    c = colorbar; ylabel(c, [params{i} ' ' units_out{i}])
    datetick('x')
    colormap(cmocean('thermal'))
    set(gca,'FontSize',15,'XLim',[min(compiled_dn) max(compiled_dn)])
  
    sp(4) = subplot(2,3,4);
    pcolor(compiled_dn(~obsNaN), -binned_depth, stats.LiveOcean.Residuals.(params{i})'); shading flat; hold on
    ylabel('Depth'); xlabel('time')
    title('Sensor - LiveOcean')
    c = colorbar; ylabel(c, ['\Delta ' units_out{i}])
    datetick('x')
    colormap(cmocean('balance'))
    set(gca,'FontSize',15,'XLim',[min(compiled_dn) max(compiled_dn)])

    sp(5) = subplot(2,3,5);
    pcolor(compiled_dn(~obsNaN), -binned_depth,stats.SSC.Residuals.(params{i})'); shading flat; hold on
    ylabel('Depth'); xlabel('time')
    title('Sensor - SalishSeaCast')
    c = colorbar; ylabel(c, ['\Delta ' units_out{i}])
    datetick('x')
    colormap(cmocean('balance'))
    set(gca,'FontSize',15,'XLim',[min(compiled_dn) max(compiled_dn)])
    
    saveas(gcf,[savedir '/figures/' siteID '_' params{i} '.png'])
end
%% ------------------------------------------------------------------------
% Format metadata & save export files
% -------------------------------------------------------------------------

fileSavePath = [savedir '/' siteID];

% Formatting for internal .nc files Global attributes
nccreate([fileSavePath '.nc'], 'time','dimensions',{'t',nT})
ncwrite([fileSavePath '.nc'], 'time',data_export_nc.time)
ncwriteatt([fileSavePath '.nc'],'time','units','s since 1970,01,01 00:00 (UTC)')

nccreate([fileSavePath '.nc'], 'depth','dimensions',{'z',nZ})
ncwrite([fileSavePath '.nc'], 'depth',data_export_nc.depth)
ncwriteatt([fileSavePath '.nc'],'depth','units','m')

ncwriteatt([fileSavePath '.nc'],'/','Site_ID',siteID);
ncwriteatt([fileSavePath '.nc'],'/','Latitude',siteLat);
ncwriteatt([fileSavePath '.nc'],'/','Longitude',siteLon);


for i = 1:numel(params)
    if strcmp(params{i}, 'time') || strcmp(params{i}, 'depth')
        continue
    end

    nccreate([fileSavePath '.nc'],[params{i} '_sensor'],'dimensions',{'t',nT,'z',nZ},'FillValue','disable')
    ncwrite([fileSavePath '.nc'],[params{i} '_sensor'],data_export_nc.([params{i} '_sensor']))
    ncwriteatt([fileSavePath '.nc'],[params{i} '_sensor'], 'units', units_out{i})
    ncwriteatt([fileSavePath '.nc'],[params{i} '_sensor'], 'source', DataURL)

    nccreate([fileSavePath '.nc'],[params{i} '_SSC'],'dimensions',{'t',nT,'z',nZ},'FillValue','disable')
    ncwrite([fileSavePath '.nc'],[params{i} '_SSC'],data_export_nc.([params{i} '_SSC']))
    ncwriteatt([fileSavePath '.nc'],[params{i} '_SSC'], 'units', units_out{i})
    ncwriteatt([fileSavePath '.nc'],[params{i} '_SSC'], 'source', sscURL)    

    nccreate([fileSavePath '.nc'],[params{i} '_LiveOcean'],'dimensions',{'t',nT,'z',nZ},'FillValue','disable')
    ncwrite([fileSavePath '.nc'],[params{i} '_LiveOcean'],data_export_nc.([params{i} '_LiveOcean']))
    ncwriteatt([fileSavePath '.nc'],[params{i} '_LiveOcean'], 'units', units_out{i})
    ncwriteatt([fileSavePath '.nc'],[params{i} '_LiveOcean'], 'source', loURL)    

    nccreate([fileSavePath '.nc'],[params{i} '_LiveOcean_interpolated'],'dimensions',{'t',nT,'z',nZ},'FillValue','disable')
    ncwrite([fileSavePath '.nc'],[params{i} '_LiveOcean_interpolated'],data_export_nc.([params{i} '_LiveOcean_interpolated']))
    ncwriteatt([fileSavePath '.nc'],[params{i} '_LiveOcean_interpolated'], 'units', units_out{i})
    ncwriteatt([fileSavePath '.nc'],[params{i} '_LiveOcean_interpolated'], 'source', 'linear interpolation of model data to match obs time stamps')  

    nccreate([fileSavePath '.nc'],[params{i} '_SSC_interpolated'],'dimensions',{'t',nT,'z',nZ},'FillValue','disable')
    ncwrite([fileSavePath '.nc'],[params{i} '_SSC_interpolated'],data_export_nc.([params{i} '_SSC_interpolated']))
    ncwriteatt([fileSavePath '.nc'],[params{i} '_SSC_interpolated'], 'units', units_out{i})
    ncwriteatt([fileSavePath '.nc'],[params{i} '_SSC_interpolated'], 'source', 'linear interpolation of model data to match obs time stamps')  
end

end