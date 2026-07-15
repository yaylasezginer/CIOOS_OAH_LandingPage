% Create a list of viable ONC sensors for incorporation into OceanECO.
% Viable sensors excludes sensors on mobile platforms or sensors deployed
% by the Community Fishers program.
%
% Uses ONC Open API (https://data.oceannetworks.ca/OpenAPI) to query sensors
% Uses https://curlconverter.com/matlab/ to convert OpenAPI curl commands
% to Matlab script

locationCodesWanted = {'PAC','ATL','ARCT','ANT'};

% Only includes sensors with propertyCodes:

propertyCodesWanted = { ...
    'co2concentration',...
    'co2partialpressure',...
    'conductivity',...
    'seawatertemperature',...
    'salinity',...
    'oxygen',...
    'ph'};


%% Web Access using Data Import and Export API

% HTTP Interface
import matlab.net.*
import matlab.net.http.*

baseURI = 'https://data.oceannetworks.ca/api/locations';

locationCode = [];
locationName = [];
propertyCode = [];
lat = [];
lon = [];
depth = [];

for l = 1:numel(locationCodesWanted)
    locCode = locationCodesWanted{l};

    for p = 1:numel(propertyCodesWanted)
        propCode = propertyCodesWanted{p};

        params = {'locationCode', locCode, 'propertyCode', propCode,'includeChildren', 'true','token', '2bc9939e-d776-41ff-abe3-3cf7b0b465d1'};

        try
            response = onc.getLocations(params);

            % Save compiled data in a new data table
            locationName = [locationName; {response.locationName}'];
            locationCode = [locationCode; {response.locationCode}'];
            propertyCode = [propertyCode; repmat({propCode},height(response),1)];
            lat = [lat; {response.lat}'];
            lon = [lon; {response.lon}'];
            depth = [depth; {response.depth}'];
            
            disp([num2str(height(response)) ' ' propCode ' detected in ' locCode])
        catch
            disp(['No ' propCode ' detected in ' locCode])
            continue
        end

    end
end

ONCsensorNetwork = table(locationCode,locationName,propertyCode,lat,lon,depth);

% Remove CF locations

CF = contains(ONCsensorNetwork.locationCode, 'CF');
ONCsensorNetwork = ONCsensorNetwork(~CF,:);

% Keep remaining unique responses

[locationKeep, ikeep, ~] = unique(ONCsensorNetwork.locationCode);
propertyCodeKeep = repmat('',height(locationKeep),1);

for i = 1:numel(locationKeep)
    ind = strcmp(ONCsensorNetwork.locationCode, locationKeep{i});
    propertyCodeKeep{i} = char(join(ONCsensorNetwork.propertyCode(ind),','));
end

ONCsensorNetwork = ONCsensorNetwork(ikeep,:);
ONCsensorNetwork.propertyCode = propertyCodeKeep';