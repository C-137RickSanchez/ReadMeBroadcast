classdef Coord_Eph 
% Coord_Eph - A class that defines the coordinates of satellites using
%             ephemerides
% Written by locateTempUserBash
% visit the user page @ github for further information
% or email using locateTempUserBash@yahoo.com
    properties
        x,y,z
        sat_clk_corr % satellite clock correction
        time@date_time        % time
        group_delay  % group delay time
        PRN          % satellite number
    end
    methods (Static)
            % getter and setter functions
            function result = get(b,par)
                     result =  b.(par);
            end
            function b = set(b,par,in)
                     b.(par)=in;
            end
           
    end    
    methods
        % The constructors are as in the following:
        % coord_obj = Coord_Eph()
        % coord_obj = Coord_Eph(x,y,z,t,PRN,sat_clk_corr,group_delay) 
        %               can be used for testing
        function coord_obj = Coord_Eph(varargin)
            na = length(varargin);
            switch na
                case 0
                    coord_obj.x=0;coord_obj.y=0;coord_obj.z=0;
                    coord_obj.PRN=0;coord_obj.sat_clk_corr=0;
                    coord_obj.time=date_time();coord_obj.group_delay=0;
                case 7
                    [xdummy, ydummy, zdummy, tdummy, sat_no, clock_corr, delay] = deal(varargin{:});
                    coord_obj.x=xdummy;coord_obj.y=ydummy;coord_obj.z=zdummy;
                    coord_obj.time=tdummy;coord_obj.PRN=sat_no;coord_obj.group_delay=delay;
                    coord_obj.sat_clk_corr=clock_corr;
                otherwise
                    na
                    error('Coord_Eph: incorrect number of arguements');
            end
          end
          function [obj, data] = broadcast_position(RINEX_fname,time_interval,coord_obj) 
            % This function reads broadcast ephemeris from a given RINEX file
            % and calculates the satellite's position for a given time interval
            %
            % Function call:
            %
            %   [obj, data] = broadcast_position(obj,RINEX_fname,time_interval)
            % 
            % Parameters and attributes:
            % 
            %   obj - ephemerides coordinate instance 
            %   data - ephemerides data 
            %   RINEX_fname - filename of the RINEX file
            %   time_interval - the time interval
            %
            % Other functions that are called:
            %
            %   broadcast_coord_all(RINEX_fname,obsTime)
            %   read_rinex(RINEX_fname)
            %
            
             % Printer that helps in debugging 
                fprintf('broadcast_position_line_68');
                fprintf('\n');
                [header,data] = read_rinex(RINEX_fname,coord_obj);  
                [n,m]=size(data);
                % Initialisation   
                dLastEpoch =[];
                dFirstEpoch = [];
                % Find the first and last Epoch 
                first_time_flag = 0;
                for i=1:n
                    for j=1:m
                        % check whether or not any data is apparent for the satellite  
                        check = arrayfun(@(x) numel(x.iEpochYear), data(i,j));
                        if check<1e-13
                           % formatSpec = 'An empty array for sat PRN:%d has been spotted.';
                           % A1 = data(i,j).iPRN;
                           % fprintf(formatSpec,A1);
                           % fprintf('\n');
                            continue;
                        end
                        BB(1) = get(data(i,j) , 'iEpochYear');
                        BB(2) = get(data(i,j) , 'iEpochMonth');                        
                        BB(3) = get(data(i,j) , 'iEpochDay');               
                        BB(4) = get(data(i,j) , 'iEpochHour');               
                        BB(5) = get(data(i,j) , 'iEpochMinute');        
                        BB(6) = get(data(i,j) , 'iEpochSecond');
                        % time to Julian date
                        dJD = date_time(BB);
                        if first_time_flag ==0
                            dFirstEpoch = dJD;
                            dLastEpoch = dJD;
                            first_time_flag=1;
                        elseif (dJD - dFirstEpoch) < 0              
                            dFirstEpoch = dJD;               
                        elseif (dJD - dLastEpoch) > 0             
                            dLastEpoch = dJD;               
                        end
                    end
                end    
            % start of the interval
                dStart = dFirstEpoch - 3*3600;
            % end of the interval
                dEnd = dLastEpoch + 3*3600;  
            %Time span in [s]
                TimeSpan = dEnd - dStart;   
            %Number of intervals
                Nint = floor(TimeSpan/time_interval);
                data(1,1) = set(data(1,1), 'dtStart', dStart);
                data(1,1) = set(data(1,1), 'dtEnd', dEnd);
            % Calculate satellite coordinates for each interval
                epo = dStart; % counter
                obj = coord_obj;
                for i = 1:Nint
                     obj.time(i) = epo;
                    % Here should the result be placed in a CoordEph class
                    vXYZ = broadcast_coord_all(data, obj);
                    [iRow,iCol] = size(vXYZ);       
                    for p=1:iRow
                            prn = p;
                            obj.PRN(prn,i) = prn;
                            obj.x(prn,i)   =  vXYZ(p,1);
                            obj.y(prn,i)   =  vXYZ(p,2);
                            obj.z(prn,i)   =  vXYZ(p,3);
                            obj.sat_clk_corr(prn,i) =  vXYZ(p,4);
                            obj.group_delay(prn,i) =  vXYZ(p,5);
                    end
                    epo = epo + time_interval;
                end
            end

            function vCoordinateVector = broadcast_coord_all(mData,coord_obj)
            % This function calculates the coordinates for all available satellites from a given RINEX file
            % 
            % Function call:
            %
            %     vCoordinateVector = broadcast_coord_all(RINEX_fname,coord_obj)
            %
            % Parameters and attributes: 
            % 
            %    mData - broadcast data pulled off from the RINEX file
            %    coord_obj  - Coord_Eph instance
            %
            % Other functions that are called:
            %   
            %   calc_coord(obs,epoch)
            %
            % 
            % DIFLIM = 3*3600; %maximum allowed time difference between toe and epoch
            
            % Printer that helps in debugging 
            fprintf('broadcast_coord_all_line_158');
            fprintf('\n');
            % Initialise  
            [n m]=size(mData);
            XYZ=[];
            vCoordinateVector = zeros(32,5);  %allocate space for 32 satellites
            epoch = coord_obj.time(end);
            dtStart = get(mData(1,1), 'dtStart');
            dtEnd = get(mData(1,1), 'dtEnd');
            
            cond1=(epoch - dtStart); 
            cond2=(epoch - dtEnd);
            
            ii=1;
            sum1=0;
            [~, s1]=size(cond1);
            while ii<=s1
                sum1=sum1+cond1(ii);
                ii=ii+1;
            end
            ii=1;
            sum2=0;
            [~, s2]=size(cond2);
            while ii<=s2
                sum2=sum2+cond2(ii);
                ii=ii+1;
            end
            if (sum1) < 0 || (sum2) > 0
                epoch
                dtStart
                dtEnd
                error('No broadcast ephemeris for this epoch');
            end
            clear ii; clear sum1; clear sum2; clear cond1; clear cond2;
            NoData = 1;
            for i=1:n  %loop over satellite
                A=[];
                MinDiff = 1e15;  %Time difference in [s]
                MinIdx = 0;  %index to the closest epoch
                for j=1:m  %loop over toe
                    % Check whether or not the satellite data exists
                    check = arrayfun(@(x) numel(x.iEpochYear), mData(i,j));
                    if abs(check) < 1e-13
                        % NoData = 1;
                        continue;
                    end
                    B(1) = get(mData(i,j) , 'iEpochYear');
                    B(2) = get(mData(i,j) , 'iEpochMonth');                        
                    B(3) = get(mData(i,j) , 'iEpochDay');               
                    B(4) = get(mData(i,j) , 'iEpochHour');               
                    B(5) = get(mData(i,j) , 'iEpochMinute');        
                    B(6) = get(mData(i,j) , 'iEpochSecond');    
                   
                    dJDSat = date_time(B); % Convert the time  
                    %dJDSat_is_date_time=isa(dJDSat,'date_time')
                    %epoch_is_a_date_time=isa(epoch,'date_time')
                    Tdif = abs(dJDSat - epoch);
                    if Tdif < MinDiff
                        MinDiff = Tdif;
                        MinIdx = j;
                        NoData = 0;
                    end
                end
            %     if MinDiff > DIFLIM  %no eph data for current epoch
            %         continue;
            %     end
                if NoData
                    continue;
                end
                [xdummy,ydummy,zdummy,ds,TGD] =  calc_coord( mData(i, MinIdx), coord_obj);
                vCoordinateVector(i,:) = [xdummy ydummy zdummy ds TGD];
                NoData = 1;
            end
             clear epoch;
            end

            function [x,y,z,ds,TGD]=calc_coord(obs,coord_obj)

            % This function calculates coordinates of one satellite at observation time
            % "epoch" from observations in the observation vector "obs". 
            % 
            % To call the function, use 
            %
            %    [x,y,z,ds,TGD]=calc_coord(obs,epoch)
            %
            % where 
            %   (x,y,z)  - are satellite coordinates in WGS84
            %   ds       - delta SV PRN code phase time offset in seconds
            %   TGD      - Differential Group delay in seconds
            %   obs     - observation vector
            %   epoch - The current observation time in reciver time
            %
            % Other functions that are called:
            %       date_time(y, mo, d, ho, min, sec) - date_time object
            %
            
            % Printer that helps in debugging 
            % fprintf('calc_coord_line_255')
            % fprintf('\n');
             epoch=coord_obj.time(end);
            % PRN/EPOCH/SC CLK
                SatNr=get(obs,'iPRN');
                EpochYear =get(obs,'iEpochYear');
                EpochMonth=get(obs,'iEpochMonth');            
                EpochDay=get(obs,'iEpochDay');            
                EpochHour=get(obs,'iEpochHour');            
                EpochMinute=get(obs,'iEpochMinute');           
                EpochSecond=get(obs,'iEpochSecond');

                dtTocl = date_time(EpochYear, EpochMonth, EpochDay, EpochHour, EpochMinute, EpochSecond);            
                a0=get(obs,'dClockBias');           
                a1=get(obs,'dClockDrift');            
                a2=get(obs,'dClockDriftRate');           
            % Broadcast item 1
                IODE = get(obs,'dIDOE');        % Issue of data (ephemeris)           
                Crs = get(obs,'dCrs');          % Amplitude of second-order harmonic pertubations           
                Delta_n = get(obs,'dDeltaN');   % Mean motion difference from computed value            
                M0 = get(obs,'dM0');            % Mean anomaly at reference time           
            % Broadcast item 2
                Cuc=get(obs,'dCuc');           % Amplitude of second-order harmonic pertubations           
                e=get(obs,'dEccent');          % Eccentricity            
                Cus=get(obs,'dCus');           % Amplitude of second-order harmonic pertubations           
                a=(get(obs,'dSqrtA'))^2;       % (Square root of the semi major axis)^2           
            % Broadcast item 3
                Toe=get(obs,'dToe');           % Ephemeris reference time           
                Cic=get(obs,'dCic');           % Amplitude of second-order harmonic pertubations            
                OMEGA=get(obs,'dOMEGA');       % Longitude of ascending node of orbit plane at beginning of week           
                Cis=get(obs,'dCis');           % Amplitude of second-order harmonic pertubations           
            % Broadcast item 4
                i0=get(obs,'di0');             % Inclination angle at reference time            
                Crc=get(obs,'dCrc');           % Amplitude of second-order harmonic pertubations           
                omega=get(obs,'dOmega');       % Argument of perigee            
                OMEGA_DOT=get(obs,'dOMEGADot');% Rate of right ascension                       
            % Broadcast item 5
                IDOT=get(obs,'dIdot');                     % Rate of incliniation angle          
                Codes_On_L2_channel=get(obs,'dCodeOnL2');            
                GPS_WEEK=get(obs,'dGpsWeek');              % To go with TOE
                dtToe = date_time(GPS_WEEK, Toe);  %convert to DateTime            
                L2_data_flag=get(obs,'dPDataFlag');        
            % Broadcast item 6
                SV_accuracy=get(obs,'dSVaccur');       % (meters)          
                SV_health=get(obs,'dSVhealth');        % (MSB only)          
                TGD=get(obs,'dTGD');                   % (seconds)   
                IODC_Issue_of_data=get(obs,'dIODC'); ; % (Clock) 
            % Broadcast item 7
                          tr_time_of_message=get(obs,'dTransTime'); % Transmission time of message (sec of GPS week))
                                                        % derived e. g. from Z-count in Hand Over Word
                                                        % (HOW)
                cons=Constants;                                        
                n0=sqrt(cons.my/a^3); % Computed mean motion - rad/sec
               
                % SV PRN code phase time offset  
                diffTcl = epoch - dtTocl;
                delta_ts = a0 + a1*diffTcl + a2*diffTcl^2; % Delta SV PRN code phase time offset in seconds
                                                            % determined the first time without relativistic 
                                                            % effects.
                % Time from ephemeris reference epoch        
                tk = epoch - dtToe - delta_ts;         
                % Corrected mean motion    
                n=n0+Delta_n;      
                % Mean anaomaly            
                Mk=M0+n*tk;           
                % Itteration to determine Kepler's equation for eccentric anomaly              
                Ek=Mk;                 
                diff=1;                   
                while abs(diff)>1.0e-13
                        diff=Mk-Ek+e*sin(Ek);           
                        Ek=Ek+diff;
                end
                % A second run up to delta_ts is initiated in order to apply the relativistic
                % effects

                % Relativistic effects 
                F=-2*sqrt(cons.my)/(cons.C)^2; 
                dTr = F * e * sqrt(a) * sin(Ek);
                % delta_ts with relativistic effects
                delta_ts=a0+a1*diffTcl + a2*diffTcl^2 + dTr; % Delta SV PRN code phase time offset in seconds
                                                            % determined the first time without relativistic 
                                                            % effects.
                % Time from ephemeris reference epoch        
                tk = epoch - dtToe - delta_ts;  
                % Itteration to determine Kepler's equation for eccentric anomaly
                % initial values

                diff=1;

                % While loop            
                while abs(diff)>1.0e-13
                        diff=Mk-Ek+e*sin(Ek);           
                        Ek=Ek+diff;
                end
                % Keppler elements
                % True anomaly            
                Cvk=(cos(Ek)-e)/(1-e*cos(Ek));           
                Svk=(sqrt(1-e^2)*sin(Ek))/(1-e*cos(Ek));            
                fk=atan2(Svk,Cvk);
                    if fk<0
                        fk=fk+2*pi;                  
                    end
                % Argument of latituide
                    fi_k=fk+omega;           
                % Second harmonic perturbations           
                    ra_uk=Cus*sin(2*fi_k)+Cuc*cos(2*fi_k);                
                    ra_rk=Crc*cos(2*fi_k)+Crs*sin(2*fi_k);                
                    ra_ik=Cic*cos(2*fi_k)+Cis*sin(2*fi_k);               
                % Corrected argument of latitude            
                    uk=fi_k+ra_uk;                        
                % Corrected radius            
                    rk=a*(1-e*cos(Ek))+ra_rk;                   
                % Corrected inclination            
                    ik=i0+ra_ik+IDOT*tk;                
                % Position in orbital plane            
                    xk=rk*cos(uk);                
                    yk=rk*sin(uk);               
                % Corrected longitude of ascending node           
                    OMEGA_k=OMEGA+(OMEGA_DOT-cons.OmegaDotE)*tk-cons.OmegaDotE*Toe;                       
                % Earth fixed coordinates
                    x=xk*cos(OMEGA_k)-yk*cos(ik)*sin(OMEGA_k);  
                    y=xk*sin(OMEGA_k)+yk*cos(ik)*cos(OMEGA_k); 
                    z=yk*sin(ik);
                    ds = delta_ts;  
                   
                 clear epoch;
            end
            
            function [A, data] = read_broadcast(coord_obj,RINEX_fname,dInterval)
            % This function reads broadcast ephemeris from a RINEX file
            % and calculates the satellite position with the wanted intervall
            %
            % To call the function use:
            %
            %      [A, data] = read_broadcast(coord_obj,sFileName,dInterval)
            %
            % where
            %   coord_obj - Coord_Eph instance 
            %   data - ephemerides data 
            %   sFileName - filename of the RINEX file
            %   dInterval - the time intervall
            %
            % Other functions that are called:
            %       read_rinex(RINEX_fname)
            
            % Printer that helps in debugging 
                fprintf('read_broadcast_line_401');
                fprintf('\n');
                [header,data] = read_rinex(RINEX_fname);
                [n m]=size(data);
            % Initialise
                dLastEpoch =[];
                dFirstEpoch = [];
            % Find the first and last Epoch 
                first_time_flag = 0;  % indicates that the process runs for the first time
                for i=1:n
                    for j=1:m
                        check = arrayfun(@(x) numel(x.iEpochYear), mData(i,j));
                        if abs(check) < 1e-13 % no data for this satellite
                            continue;
                        end
                        BB(1) = get(data(i,j) , 'iEpochYear');
                        BB(2) = get(data(i,j) , 'iEpochMonth');                        
                        BB(3) = get(data(i,j) , 'iEpochDay');               
                        BB(4) = get(data(i,j) , 'iEpochHour');               
                        BB(5) = get(data(i,j) , 'iEpochMinute');        
                        BB(6) = get(data(i,j) , 'iEpochSecond');    
                        %if abs(BB(1)) < 1e-13 % no data for this satellite
                        %    continue;
                        %end
                        dJD = DateTime(BB); % Convert the time to JD
                        %if i==1 & j==1  %ORIGINAL
                        if first_time_flag ==0
                            dFirstEpoch = dJD;
                            dLastEpoch = dJD;
                            first_time_flag=1;
                        elseif (dJD - dFirstEpoch) < 0  
                            dFirstEpoch = dJD;
                        elseif (dJD - dLastEpoch) > 0  
                            dLastEpoch = dJD;
                        end  
                    end
                end    
            % Define start of the interval
                dStart = dFirstEpoch - 3*3600;
            % Define the end of the interval
                dEnd = dLastEpoch + 3*3600;
            %Time span in [s]
                TimeSpan = dEnd - dStart;   
            %Number of intervals    
                Nint = floor(TimeSpan/dInterval);    
                data(1,1) = set(data(1,1), 'dtStart', dStart);
                data(1,1) = set(data(1,1), 'dtEnd', dEnd);
            % Calculate satellite coordinates for each interval
                epo = dStart; % counter
                for i = 1:Nint
                    coord_obj.dTime(i) = epo;
                    % Here should the resulult be placed in a CoordEph class
                    vXYZ = beEphCoordCalc(data, epo);
                    [iRow,iCol] = size(vXYZ);       
                    for p=1:iRow
                            prn = p;
                            coord_obj.iPRN(prn,i) = prn;
                            coord_obj.x(prn,i)   =  vXYZ(p,1);
                            coord_obj.y(prn,i)   =  vXYZ(p,2);
                            coord_obj.z(prn,i)   =  vXYZ(p,3);
                            coord_obj.sat_clk_corr(prn,i) =  vXYZ(p,4);
                            coord_obj.group_delay(prn,i) =  vXYZ(p,5);
                    end
                    epo = epo +  dInterval;
                end
            end
            function [Header,out] = read_rinex(RINEX_fname,coord_obj)
            % This function reads the header and ephemeris from RINEX format
            %
            % The function gives/takes the following output/input parameters 
            %
            %    [Header,out] = read_rinex(RINEX_fname)
            %
            % where
            %    RINEX_fname - RINEX file name
            %    Header  - Return the header of the RINEX file
            %    out - Returns the data from the RINEX file
            % 
            % Other functions that are called:
            %   RINEX_header - ephemerides header  
            %
            
            % Printer that helps in debugging
                fprintf('read_rinex_line_484');
                fprintf('\n');
            % Opnening the RINEX file
                fid=fopen(RINEX_fname);
            % Initialise
                rad=0; % Start value Counter
                num=0; % Start value Counter
                tline= fgetl(fid);
                countComments=0;
                sat_list=[]; %  a vector with all the satellites in the file (see below)
                A=RINEX_header; % Creates a Header object
                while num<1 % Loop until END OF HEADER
                    matches = findstr(tline,'END OF HEADER');
                    num = length(matches);
                    if num==0; % If header is not on the line take next
                        % RINEX VERSION /Type of observation
                        matches = findstr(tline,'RINEX VERSION');
                        num = length(matches);
                        if num>0
                            A = set(A,'iRinexVersion',str2double(tline(6:(6+14))));
                            A = set(A,'sFileType',tline(21:40));
                            num=0;
                        end
                        % PGM / RUN BY / DATE
                        matches = findstr(tline,'PGM / RUN BY / DATE');
                        num = length(matches);
                        if num>0
                            A =set(A,'sPGM',tline(1:20));
                            A =set(A,'sRunBy',tline(21:40));
                            A =set(A,'sDate',tline(41:60));
                            num=0;    
                        end
                        % COMMENTS Comment line(s)  
                        matches = findstr(tline,'COMMENT');
                        num = length(matches);
                        if num>0
                            countComments=countComments+1;
                            if countComments==1;
                                sComment=tline(1:60);
                            else 
                                sComment=[A.aComment;tline(1:60)];
                            end
                            A =set(A,'sComment',sComment);
                            num=0;
                        end
                        % ION ALPHA, Ionosphere parameters A0-A3 of almanac
                        matches = findstr(tline,'ION ALPHA');
                        num = length(matches);
                        if num>0
                            A=set(A,'dAlfaIon1',str2double(strrep(tline( 3:15),'D','e')));
                            A=set(A,'dAlfaIon2',str2double(strrep(tline(16:28),'D','e')));
                            A=set(A,'dAlfaIon3',str2double(strrep(tline(29:40),'D','e')));
                            A=set(A,'dAlfaIon4',str2double(strrep(tline(41:54),'D','e')));
                            num=0;    
                        end
                        % ION BETA, Ionosphere parameters B0-B3 of almanac
                        matches = findstr(tline,'ION BETA');
                        num = length(matches);
                        if num>0
                            A=set(A,'dBetaIon1',str2double(strrep(tline( 3:15),'D','e')));
                            A=set(A,'dBetaIon2',str2double(strrep(tline(16:28),'D','e')));
                            A=set(A,'dBetaIon3',str2double(strrep(tline(29:40),'D','e')));
                            A=set(A,'dBetaIon4',str2double(strrep(tline(41:54),'D','e')));
                            num=0;    
                        end
                        % DELTA-UTC, Almanac parameters to compute time in UTC
                        matches = findstr(tline,'DELTA-UTC');
                        num = length(matches);
                        if num>0
                           A=set(A,'dA0',str2double(strrep(tline( 4:23),'D','e')));
                           A=set(A,'dA1',str2double(strrep(tline(24:43),'D','e')));
                           A=set(A,'dRefTime',str2double(strrep(tline(44:53),'D','e')));
                           A=set(A,'dRefWeek',str2double(strrep(tline(54:60),'D','e')));
                           num=0;    
                        end
                        % LEAP SECONDS, Delta time due to leap seconds   
                        matches = findstr(tline,'LEAP SECONDS');
                        num = length(matches);
                        if num>0
                            A=set(A,'iLeapSeconds',str2double(strrep(tline(1:6),'D','e')));
                            num=0; 
                        end
                    tline= fgetl(fid);
                    rad=rad+1; % Counter
                    end
                end
                Header = A;
                C = ephemeris_data;
                while feof(fid) == 0; % Loop continues until the end of the file
                     B=ephemeris_data;
                    % PRN / EPOCH / SV CLK
                        tline= fgetl(fid); rad=rad+1; % Take a new line from the file
                        B=set(B,'iPRN',str2double(strrep(tline(1:2),'D','e')));
                        year=str2double(tline(3:5));
                            % Convert the year to 4 digits
                            if year<94;
                                year=year+2000;
                            else
                                year=year+1900;
                            end 
                        B=set(B,'iEpochYear',year);
                        B=set(B,'iEpochMonth',str2double(strrep(tline(6:8),'D','e')));    
                        B=set(B,'iEpochDay',str2double(strrep(tline(9:11),'D','e')));
                        B=set(B,'iEpochHour',str2double(strrep(tline(12:14),'D','e')));
                        B=set(B,'iEpochMinute',str2double(strrep(tline(15:18),'D','e')));        
                        B=set(B,'iEpochSecond',str2double(strrep(tline(19:22),'D','e')));    
                        B=set(B,'dClockBias',str2double(strrep(tline(23:41),'D','e')));
                        B=set(B,'dClockDrift',str2double(strrep(tline(42:60),'D','e')));
                        B=set(B,'dClockDriftRate',str2double(strrep(tline(61:79),'D','e')));   
                    % Checking if there is any observation data in the instance
                    % and return with a boolean in order to inform that
                    % there actually is data inside the instance
                    % Broadcast line 1
                        tline= fgetl(fid);rad=rad+1;
                        B=set(B,'dIDOE',str2double(strrep(tline(4:22),'D','e')));    
                        B=set(B,'dCrs',str2double(strrep(tline(23:41),'D','e')));
                        B=set(B,'dDeltaN',str2double(strrep(tline(42:60),'D','e')));
                        B=set(B,'dM0',str2double(strrep(tline(61:79),'D','e')));  
                    % Broadcast line 2
                        tline= fgetl(fid);rad=rad+1;
                        B=set(B,'dCuc',str2double(strrep(tline(4:22),'D','e')));    
                        B=set(B,'dEccent',str2double(strrep(tline(23:41),'D','e')));
                        B=set(B,'dCus',str2double(strrep(tline(42:60),'D','e')));
                        B=set(B,'dSqrtA',str2double(strrep(tline(61:79),'D','e')));
                    % Broadcast line 3
                        tline= fgetl(fid);rad=rad+1;
                        B=set(B,'dToe',str2double(strrep(tline(4:22),'D','e')));    
                        B=set(B,'dCic',str2double(strrep(tline(23:41),'D','e')));
                        B=set(B,'dOMEGA',str2double(strrep(tline(42:60),'D','e')));
                        B=set(B,'dCis',str2double(strrep(tline(61:79),'D','e')));
                     % Broadcast line 4
                        tline= fgetl(fid);rad=rad+1;
                        B=set(B,'di0',str2double(strrep(tline(4:22),'D','e')));    
                        B=set(B,'dCrc',str2double(strrep(tline(23:41),'D','e')));
                        B=set(B,'dOmega',str2double(strrep(tline(42:60),'D','e')));
                        B=set(B,'dOMEGADot',str2double(strrep(tline(61:79),'D','e')));
                    % Broadcast line 5
                        tline= fgetl(fid);rad=rad+1;
                        B=set(B,'dIdot',str2double(strrep(tline(4:22),'D','e')));    
                        B=set(B,'dCodeOnL2',str2double(strrep(tline(23:41),'D','e')));
                        B=set(B,'dGpsWeek',str2double(strrep(tline(42:60),'D','e')));
                        B=set(B,'dPDataFlag',str2double(strrep(tline(61:79),'D','e')));
                    % Broadcast line 6
                        tline= fgetl(fid);rad=rad+1;
                        B=set(B,'dSVaccur',str2double(strrep(tline(4:22),'D','e')));    
                        B=set(B,'dSVhealth',str2double(strrep(tline(23:41),'D','e')));
                        B=set(B,'dTGD',str2double(strrep(tline(42:60),'D','e')));  
                        B=set(B,'dIODC',str2double(strrep(tline(61:79),'D','e')));
                    % Broadcast line 7
                        tline= fgetl(fid);rad=rad+1;
                        if length(tline)>78
                            B=set(B,'dTransTime',str2double(strrep(tline(4:22),'D','e')));
                            B=set(B,'dSpare1',str2double(strrep(tline(23:41),'D','e')));
                            B=set(B,'dSpare2',str2double(strrep(tline(42:60),'D','e')));
                            B=set(B,'dSpare3',str2double(strrep(tline(61:79),'D','e')));
                        elseif length(tline)>55
                            B=set(B,'dTransTime',str2double(strrep(tline(4:22),'D','e')));
                            B=set(B,'dSpare1',str2double(strrep(tline(23:41),'D','e')));
                            B=set(B,'dSpare2',str2double(strrep(tline(42:60),'D','e')));
                        elseif length(tline)>23
                            B=set(B,'dTransTime',str2double(strrep(tline(4:22),'D','e')));
                            B=set(B,'dSpare1',str2double(strrep(tline(23:41),'D','e')));
                        elseif length(tline)<23
                            B=set(B,'dTransTime',str2double(strrep(tline(4:22),'D','e')));
                        end
                     % The following routine prevents data over-write when data from
                     % a satellite is received twice
                     % Initialise
                       state=0;       % switch
                        i=0;        % Reset counter
                     % Actual satellite
                        xx=get(B,'iPRN');
                       if isempty(sat_list) % First satellite
                            sat_list=xx;            % z - new vector with the satellite numbers
                            C(xx,1) = B;       
                            zCount=1;
                        else % otherwise
                            while ((state==0) & (i<length(sat_list))) 
                            i=i+1; % Counter
                              % check if the satellite data already is imported 
                              if xx==sat_list(i)
                                  state=1;
                              else
                                  state=0;
                              end
                            end
                            % switch case - if data not exists, enter it on the first position
                            % otherwise in the second position
                            switch state
                                case 0
                                    sat_list=[sat_list,xx];
                                    zCount=[zCount,1];
                                    C(xx,1)=B;
                                case 1
                                    zCount(i)=zCount(i)+1;
                                    C(xx,zCount(i))=B;                  
                            end            
                        end          
                        clear B;
                end
                % Return the the header and all data
                out = C;
            end
    end
end
      

function header = RINEX_header(varargin)
        % Navigation file header data, stored in the following parameters
        %
        % RINEX VERSION         : iRinexVersion
        % File Type             : sFileType
        % Name of Program       : sPGM
        % Name of Agency        : sRunBy
        % Date of File Creation : sDate
        % Comment               : sComment
        % Ionosphere param. A0  : dAlfaIon1
        % Ionosphere param. A1  : dAlfaIon2
        % Ionosphere param. A2  : dAlfaIon3
        % Ionosphere param. A3  : dAlfaIon4
        % Ionosphere param. B0  : dBetaIon1
        % Ionosphere param. B1  : dBetaIon2
        % Ionosphere param. B2  : dBetaIon3
        % Ionosphere param. B3  : dBetaIon4
        % Delta-UTC A0          : dA0
        % Delta-UTC A1          : dA1
        % Reference Time        : dRefTime
        % Reference Week        : dRefWeek
        % Leap second           : iLeapSeconds
        %
        fprintf('rinex_header_line713');
        fprintf('\n');
        na=length(varargin);
            switch na
                case 0
                    % no arguments
                    header.iRinexVersion=0;header.sFileType=0;header.sPGM=0;header.sRunBy=0;
                    header.sDate=0;header.sComment=0;header.dAlfaIon1=0;header.dAlfaIon2=0;header.dAlfaIon3=0;
                    header.dAlfaIon4=0;header.dBetaIon1=0;header.dBetaIon2=0;header.dBetaIon3=0;header.dBetaIon4=0;
                    header.dA0=0;header.dA1=0;header.dRefTime=0;header.dRefWeek=0;header.iLeapSeconds=0;
                case 1
                    header.iRinexVersion=data_in.iRinexVersion;
                    header.sFileType=data_in.sFileType;
                    header.sPGM=data_in.sPGM;
                    header.sRunBy=data_in.sRunBy;
                    header.sDate=data_in.sDate;
                    header.sComment=data_in.sComment;
                    header.dAlfaIon1=data_in.dAlfaIon1;
                    header.dAlfaIon2=data_in.dAlfaIon2;
                    header.dAlfaIon3=data_in.dAlfaIon3;
                    header.dAlfaIon4=data_in.dAlfaIon4;
                    header.dBetaIon1=data_in.dBetaIon1;
                    header.dBetaIon2=data_in.dBetaIon2;
                    header.dBetaIon3=data_in.dBetaIon3;
                    header.dBetaIon4=data_in.dBetaIon4;
                    header.dA0=data_in.dA0;
                    header.dA1=data_in.dA1;
                    header.dRefTime=data_in.dRefTime;
                    header.dRefWeek=data_in.dRefWeek;
                    header.iLeapSeconds=data_in.iLeapSeconds;
            end
end
function ephemeris = ephemeris_data()

% Broadcast Emepheris from Rinex file 
% PRN / EPOCH / SV CLK
% Satellite PRN number           : iPRN
% Epoch year                     : iEpochYear
% Epoch month                    : iEpochMonth
% Epoch day                      : iEpochDay
% Epoch hour                     : iEpochHour
% Epoch minute                   : iEpochMinute
% Epoch second                   : iEpochSecond
% Sat Clock Bias (sec)           : dClockBias
% Sat Clock drift(sec/sec)       : dClockDrift
% Sat Clock drift rate (sec/sec2): dClockDriftRate
%
% BROADCAST ORBIT - line 1
% IDOE Issue of Data, Ephemeris  : dIDOE
% Crs (meters)                   : dCrs
% Delta n (radians/sec)          : dDeltaN
% M0 (radians)                   : dM0
%
% BROADCAST ORBIT - line 2
% Cuc (radians)                  : dCuc
% e Eccenricity                  : dEccent
% Cus (radians)                  : dCus
% sqrt(A) (sqrt(m))              : dSqrtA
%
% BROADCAST ORBIT - line 3
% Toe Time of Ephemeris (sec of GPS week): dToe
% Cic (radians)                  : dCic
% OMEGA (radians)                : dOMEGA
% Cis (radians)                  : dCis
%
% BROADCAST ORBIT - line 4
% i0                             : di0
% Crc (radians)                  : dCrc
% omega (radians)                : dOmega
% OMEGA Dot (radians)            : dOMEGADot
%
% BROADCAST ORBIT - line 5
% Idot                           : dIdot
% Codes on L2 channel            : dCodeOnL2
% GPS Week # (to go with TOE)    : dGpsWeek
% L2 P data flag                 : dPDataFlag
%
% BROADCAST ORBIT - line 6
% SV Accuracy (meters)           : dSVaccur
% SV health   (MSB only)         : dSVhealth
% TGD                            : dTGD
% IODC Issue of Data, Clock      : dIODC
%
% BROADCAST ORBIT - line 7
% Transmission time of message   : dTransTime
% Spare1                         : dSpare1
% Spare2                         : dSpare2
% Spare3                         : dSpare3

ep = date_time;
ephemeris.iPRN=0;

           ephemeris.iEpochYear=0;
           ephemeris.iEpochMonth=0;
           ephemeris.iEpochDay=0;
           ephemeris.iEpochHour=0;
           ephemeris.iEpochMinute=0;
           ephemeris.iEpochSecond=0;
           ephemeris.dClockBias=0;
           ephemeris.dClockDrift=0;
           ephemeris.dClockDriftRate=0;
           ephemeris.dIDOE=0;
           ephemeris.dCrs=0;
           ephemeris.dDeltaN=0;
           ephemeris.dM0=0;
           ephemeris.dCuc=0;
           ephemeris.dEccent=0;
           ephemeris.dCus=0;
           ephemeris.dSqrtA=0;
           ephemeris.dToe=0;
           ephemeris.dCic=0;
           ephemeris.dOMEGA=0;
           ephemeris.dCis=0;
           ephemeris.di0=0;
           ephemeris.dCrc=0;
           ephemeris.dOmega=0;
           ephemeris.dOMEGADot=0;
           ephemeris.dIdot=0;
           ephemeris.dCodeOnL2=0;
           ephemeris.dGpsWeek=0;
           ephemeris.dPDataFlag=0;
           ephemeris.dSVaccur=0;
           ephemeris.dSVhealth=0;
           ephemeris.dTGD=0;
           ephemeris.dIODC=0;
           ephemeris.dTransTime=0;
           ephemeris.dSpare1=0;
           ephemeris.dSpare2=0;
           ephemeris.dSpare3=0;
           ephemeris.dtStart=ep;
           ephemeris.dtEnd=ep;
   

end
 function result = get(b,par)
      result =  b.(par);
  end
  function b = set(b,par,in)
         b.(par)=in;
  end
         