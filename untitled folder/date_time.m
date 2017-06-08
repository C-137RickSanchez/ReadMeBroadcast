
classdef date_time
% date_time - A class that defines the date and performs conversions between
%            different forms
% Written by locateTempUserBash
% visit the user page @ github for further information
% or email using locateTempUserBash@yahoo.com
    properties 
         year, month, day, hour, min, sec % date and time
         gweek % GPS week
         wsec % seconds of GPS week 
         MJD % modified julian date
         dweek % day of GPS week
         DOY % day of year
    end

   methods (Static)
       % getter and setter functions     
       function result = get(b,par)
                result = b.(par);
        end

        function result = set(b,par,in)
                 b.(par) = in;
                 result = b.(par);
        end

        function disp(A)

            fprintf('Year                     : %4d\n', A.year);
            fprintf('Month                    : %4d\n', A.month);
            fprintf('Day                      : %4d\n', A.day);
            fprintf('Hour                     : %4d\n', A.hour);
            fprintf('Minute                   : %4d\n', A.min);
            fprintf('Seconds                  : %.5f\n', A.sec);
            fprintf('GPS week                 : %4d\n', A.gweek);
            fprintf('Seconds of GPS week      : %.5f\n', A.wsec);
            fprintf('Day of week              : %4d\n', A.dweek);
            fprintf('MJD                      : %.5f\n', A.MJD);
            fprintf('Day of Year              : %.5f\n', A.DOY);
        end
   end

   methods
        % The constructors are as in the following
        % A = date_time()
        % A = date_time(MJD) 
        % A = date_time(gw, ws) 
        % A = date_time(y, mo, d, ho, min, sec)
        % A = date_time([y mo d ho min sec])
        function A = date_time(varargin)
        % This constructor includes the following functions:
        %       MJD2Date(A) - converts julian date to datum
        %       GPS2Date(A) - converts GPS time to datum
        %       Date2GPS(A) - converts datum to GPS time
         na = length(varargin);
            switch na
                case 0
                    A.MJD=0;A.gweek=0;A.dweek=0;A.wsec=0;A.year=0;A.month=0;A.day=0;A.min=0;A.sec=0;A.DOY=0;
                case 6  %input in georgian datum and time of day
                    [y, mo, d, h, m, s] = deal(varargin{:});
                    A.MJD=0;A.gweek=0;A.dweek=0;A.wsec=0;
                    A.year=y;
                    A.month=mo;
                    A.day=d;
                    A.hour=h;
                    A.min=m;
                    A.sec=s;
                    A.DOY=0;
                    A = date2GPS(A);
                case 2  %input in GPS week and wsec
                    [gw, ws] = deal(varargin{:});
                    A.gweek=gw;A.wsec=ws;
                    A.MJD=0;A.dweek=0;A.year=0;A.month=0;A.day=0;A.hour=0;A.min=0;A.sec=0;A.DOY=0;
                    A = GPS2date(A);
                case 1 
                    n = length(varargin{1});
                    if n == 1 %input in MJD
                        mjd = varargin{1};
                        A.MJD=mjd;A.gweek=0;A.dweek=0;A.wsec=0;A.year=0;A.month=0;A.day=0;A.hour=0;A.min=0;A.sec=0;A.DOY=0;
                        A = MJD2date(A);
                    elseif n >= 6
                        in=varargin{1};
                        A.year = in(1);
                        A.month = in(2);
                        A.day = in(3);
                        A.hour = in(4);
                        A.min = in(5);
                        A.sec = in(6);
                        A = date2GPS(A);
                    end
                    
                otherwise
                    na
                    error('date_time: Incorrect number of arguments');
            end
        end
        function  A = date2GPS(A)
            %date2GPS converts datum to GPS time and MJD
            %dtA = date2GPS(dtA)
            %dtA - object of type date_time
            y = A.year; 
            mo = A.month;
            if mo <= 2 
                y = y - 1; 
                mo = mo + 12;
            end
            a = 365.25*y;
            b = (mo+1)*30.6001;
            dh = A.hour + A.min/60 + A.sec/3600;  %hours in day
            jd = floor(a) + floor(b) + A.day + 1720981.5;  %+ dh/24  
            A.MJD = jd-2400000.5 + dh/24; 
            a = (jd - 2444244.5)/7;
            A.gweek = floor(a);
            wsec = (a - A.gweek)*7.*86400.;         % seconds of the week - not sufficient precision
            A.dweek = round(wsec/86400.);
            A.wsec = A.dweek*86400 + dh*3600;     % seconds of the week -  sufficient precision
        end
        function A = GPS2date(A)
            %GPS2date converts datum to GPS time and MJD
            %A = GPS2date(A)
            %A - object of type date_time
            jd = A.gweek*7 + A.wsec/86400 + 2444244.5;
            A.MJD = jd-2400000.5;
            a = floor(jd+0.5); 
            b = a + 1537; 
            c = floor((b-122.1)/365.25);
            d = floor(365.25*c); 
            e = floor((b-d)/30.6001);
            f = jd+0.5;
            A.day = b - d - floor(30.6001*e); % + (int) modf(jd+0.5,&pom);
            A.month = e-1-12* floor(e/14);
            A.year = c - 4715 - floor((7+ A.month)/10);
            A.dweek = floor(A.wsec/86400.);
            pom = A.wsec/3600 - A.dweek*24;
            A.hour = floor(pom);
            pom = (pom - A.hour)*60;
            % A.min = floor(pom);
            % A.sec = (pom - A.min)*60;
            A.sec = A.wsec - A.dweek*86400 - A.hour*3600;  % - A.min*60;
            A.min = floor(A.sec/60);
            A.sec = A.sec - A.min*60;           
        end
        function [doy, ep] = dayofyear(ep)

            %computes day of year

            %y = get(ep,'year');
            ep0 = date_time(ep.year,1,1,0,0,0);
            %doy = floor(get(ep,'MJD') - get(ep0,'MJD'))+1;
            ep.DOY = floor(ep.MJD - ep0.MJD) + 1;
            doy = ep.DOY;
        end

        function r = minus(a,b)
            % Computes differences between two date_time 
            % instances of a and b
            % or instance a minus b seconds
            % the results is in seconds 
            r=[];
            if isa(a, 'date_time') & isa(b, 'date_time')
                [k,l] = size(a);
                [m,n] = size(b);
                if (k==m & l==n)
                    dw = a.gweek - b.gweek;
                    ds = a.wsec - b.wsec;
                    r = dw*86400*7 + ds;
                elseif (k==1 & l>1) & (m ==1 & n==1)
                    for j=1:l
                        dw(j) = a(j).gweek - b.gweek;
                        ds(j) = a(j).wsec - b.wsec;
                        r(j) = dw(j)*86400*7 + ds(j);
                    end
                end                           
            elseif isa(a, 'date_time') & isa(b, 'double')
                pom = a.wsec - b;
                r = date_time(a.gweek, pom);
            else
                err = sprintf('Operator minus in date_time does not allow arguments %s %s', class(a), class(b));
                error(err);
            end
            
        end

        function r = plus(a,b)
            % Computes the adition of two date_time
            % instances a and b
            % or instance a plus b seconds
            % the results is in seconds 

            if isa(a, 'date_time') && isa(b, 'double')
                pom = a.wsec + b;
                r = date_time(a.gweek, pom);
            elseif isa(b, 'date_time') && isa(a, 'double')
                pom = a + b.wsec;
                r = date_time(b.gweek, pom);
            else
                err = sprintf('Operator plus in date_time does not allow arguments %s %s', class(a), class(b));
                error(err);
            end
            
        end
        
        function A = MJD2date(A)
            %MJD2date converts MJD to GPS time
            %A = MJD2date(A)
            %A is date_time object

            jd = A.MJD + 2400000.5;
            week = (jd - 2444244.5)/7;
            A.gweek = floor(week);
            A.wsec = (week - A.gweek)*86400*7;
            A = GPS2date(A);   
        end
        
           
   end
end