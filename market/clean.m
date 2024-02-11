clc
clear
close all
load pendatanomiss.mat
Datasize=height(Austin_data);
uni_userid=unique(Austin_data.dataid);
uni_usernum=length(uni_userid);

dataid=Austin_data.dataid;
time=Austin_data.TIME;
grid=Austin_data.grid;
solar=Austin_data.solar;
for i=1:uni_usernum  % find index
   user_index=find(dataid==uni_userid(i));
   time_org=time(user_index);
   %% sort time
   [time_med,sort_index]=sort(time_org);
   dataid_med=dataid(user_index);
   grid_med=grid(user_index);
   solar_med=solar(user_index);
   combine=[dataid_med,grid_med,solar_med];
   com_med=combine(sort_index,:);
   eval(['ID_',num2str(uni_userid(i)),'=com_med;']);
   eval(['TimeID_',num2str(uni_userid(i)),'=[time_med];']);
end
for i=1:uni_usernum 
% for i=11:11
    i
   eval(['idlength= length(ID_',num2str(uni_userid(i)),')']);
   eval([' timeseries=(TimeID_',num2str(uni_userid(i)),')']);
   eval([' dataseries=(ID_',num2str(uni_userid(i)),')']);
   newdata=dataseries;
   newtime=timeseries;
   nstepcount=0;
    for j=1:idlength-1
        
        timegap=datetime(timeseries(j,1))-datetime(timeseries(j+1,1));
        if minutes(timegap)== -15
            continue;
        else 
          nstep=minutes(timegap)/(-15)-1;
          if j>nstep    
          % normal situation
          %
          datanstep_front=dataseries(j-nstep+1:j,:);
          datanstep_rear=dataseries(j+nstep+1:j+2*nstep,:); 
          timenstep_front=timeseries(j-nstep+1:j,:);
          timenstep_rear=timeseries(j+nstep+1:j+2*nstep,:); 
          datainsert=(datanstep_front+datanstep_rear)/2;
          dateser=minutes(15:15:-minutes(timegap)-15);
          timeinsert=repmat(datetime(timeseries(j,1)),1,length(dateser))+dateser;
          
          %% insert
          newdata_upper=newdata(1:j+nstepcount,:);
          newdata_down=newdata(j+nstepcount+1:end,:);
          newtime_upper=newtime(1:j+nstepcount,:);
          newtime_down=newtime(j+nstepcount+1:end,:);
          newdata=[newdata_upper;datainsert;newdata_down];
          newtime=[newtime_upper;timeinsert';newtime_down];
          nstepcount=nstep+nstepcount;
          % normal sistuation ended
          else  %filling using latter data directly
          datanstep_rear=dataseries(j+nstep+1:j+2*nstep,:); 
          timenstep_rear=timeseries(j+nstep+1:j+2*nstep,:); 
          datainsert=datanstep_rear;
          dateser=minutes(15:15:-minutes(timegap)-15);
          timeinsert=repmat(datetime(timeseries(j,1)),1,length(dateser))+dateser;
            %% insert
          newdata_upper=newdata(1:j+nstepcount,:);
          newdata_down=newdata(j+nstepcount+1:end,:);
          newtime_upper=newtime(1:j+nstepcount,:);
          newtime_down=newtime(j+nstepcount+1:end,:);
          newdata=[newdata_upper;datainsert;newdata_down];
          newtime=[newtime_upper;timeinsert';newtime_down];
          nstepcount=nstep+nstepcount;  
              
          end 
        end
%          if j==(idlength-1)&& timeseries(end,1)~=TimeID_2818(end,1)
%           timegap=datetime(timeseries(end,1))-datetime(TimeID_2818(end,1));
%           nstep=minutes(timegap)/(-15)-1;
%           %% move the previous day data
%           datanstep_rear=dataseries(j+nstep+1:j+2*nstep,:); 
%           timenstep_rear=timeseries(j+nstep+1:j+2*nstep,:); 
%           datainsert=datanstep_rear;
%           dateser=minutes(15:15:-minutes(timegap)-15);
%           timeinsert=repmat(datetime(timeseries(j,1)),1,length(dateser))+dateser;
%             %% insert
%           newdata_upper=newdata(1:j+nstepcount,:);
% 
%           newtime_upper=newtime(1:j+nstepcount,:);
%           newtime_down=newtime(j+nstepcount+1:end,:);
%           newdata=[newdata_upper;datainsert;newdata_down];
%           newtime=[newtime_upper;timeinsert';newtime_down];
%           nstepcount=nstep+nstepcount;  
%         
%              
%          end
    end
   eval([' NewTimeID_',num2str(uni_userid(i)),'=newtime;']);
   eval([' NewID_',num2str(uni_userid(i)),'=newdata;']);
end