function [pos, vel] = readSP3(filename)
%
% [pos vel] = readSP3(filename)
%
%REQUIRES
% Requires the name and path of an .sp3 file, either in full format or
% without the header as output by the SPPLEO program from DLR. This program
% ignores all records for non-GPS satellites (G01, G02 are fine but E01,
% R01, L01... are ignored) and assumes the SP3 file is using the GPS
% timescale for the time stamps.
%PROMISES
% Outputs two matrices:
% pos is an m x 7 matrix of every position record found in the file.
%   Columns are PRN, GPSweek,GPSsecond, posx, posy, posz, clock. 
%   Units are assumed to be km for position and microseconds for clock
%   offset based on the .sp3 file format definition. 
% vel is an m x 7 matrix of every velocity record found in the file.
%   Columns are PRN, GPSweek, GPSsecond, velx, vely, velz, clockdrift. 
%   Units are assumed to be dm/s for velocity and 10^-4 microseconds/second
%   for clock drift based on the .sp3 file format definition. If no
%   velocity records are present, the vel matrix is emtpy.
% This program ignores all lines which are not *, P or V records, or EOF,
% and ignores any covariance or other information in the SP3 file.
%
% Erin Kahr, Jan 2011

%initialize
filename;
fid = fopen(filename);

date = [];
pos = [];
vel = [];

if (fid == -1)
    fprintf('File %s could not be read, continuing without it!',filename);
    return;
end

temp = textscan(fid,'%s','Delimiter','\n','MultipleDelimsAsOne',1);
fclose(fid);

dateindex = 0;
posindex = 0;
velindex = 0;
for n = 1:(size(temp{1},1))

    line = temp{1}{n};
    % fprintf([line '\n']);
    
    if(strcmp(line(1),'*'))
        
        dateindex = dateindex + 1;
        date = [date;nan(1,6)];
        
        date(dateindex,1) = str2double(line(4:7));%year
        date(dateindex,2) = str2double(line(9:10));%month
        date(dateindex,3) = str2double(line(12:13));%day
        date(dateindex,4) = str2double(line(15:16));%hour
        date(dateindex,5) = str2double(line(18:19));%minute
        date(dateindex,6) = str2double(line(21:31));%second
        
        tempdate = datenum(date(dateindex,:));
        [GPSweek,GPSsecond] = datenum2GPS(tempdate,0);%UTC offset is zero because sp3 files are already in the GPS timescale
        GPSsecond = round(GPSsecond);

    elseif(strcmp(line(1),'P'))
        
       if(~strcmp(line(2),'G'))
           continue;
       end
        
       posindex = posindex + 1;
       pos = [pos; nan(1,7)];
       
       pos(posindex,1) = str2double(line(3:4));%satid/prn, ideally should be a char with fields 2:4
       pos(posindex,2) = GPSweek;
       pos(posindex,3) = GPSsecond;
       pos(posindex,4) = str2double(line(5:18));%x coordinate (km)
       pos(posindex,5) = str2double(line(19:32));%y coordinate (km)
       pos(posindex,6) = str2double(line(33:46));%z coordinate (km)
       pos(posindex,7) = str2double(line(47:60));%clock (microseconds)
       
    elseif(strcmp(line(1),'V'))
        
       if(~strcmp(line(2),'G'))
           continue;
       end
        
       velindex = velindex + 1;
       vel = [vel; nan(1,7)];
       
       vel(velindex,1) = str2double(line(3:4));%satid/prn, ideally should be a char with fields 2:4
       vel(velindex,2) = GPSweek;
       vel(velindex,3) = GPSsecond;
       vel(velindex,4) = str2double(line(5:18));%x coordinate (dm/s)
       vel(velindex,5) = str2double(line(19:32));%y coordinate (dm/s)
       vel(velindex,6) = str2double(line(33:46));%z coordinate (dm/s)
       vel(velindex,7) = str2double(line(47:60));%clock drift(10^-4 microseconds/second)

    elseif(strcmp(line(1),'E'))
        if(strcmp(line(1:3),'EOF'))
            break;
            %otherwise covariance information is ignored at the moment
        end
    end

end

function [week, second] = datenum2GPS(datenumin, utcoffset)
%
% [week second] = datenum2GPS(datenumin, utcoffset)
%
% Converts at datenum, assumed to be in UTC, to GPS week and second
%
% Erin Kahr, July 22, 2014

datenumin = datenumin+utcoffset;%GPS timescale

GPSdawn = datenum([1980 01 06 0 0 0]);

continuousdays = datenumin-GPSdawn;

week = floor(continuousdays/7);

second = continuousdays - week*7;%fractional week in days
second = second*86400;%fractional week in seconds
