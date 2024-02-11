%BTMdis
function [capacity MAPE L_dis stdpv]=btmdis(PV,NTL)
     onepv=PV(1:8760,:)';
     aggpv=sum(onepv,1);
     stdpv=aggpv/max(aggpv);
%     onedata=NTL(1:6720,:)';
     %stdpv=(sum(onepv,1)/10)/max((sum(onepv,1)/10));
     PV_std=reshape(stdpv,24,365);
     K_means= kmeans(PV_std',4);[row, col] = find(isnan(K_means));
        l=1;m=1;n=1;o=1;
        for i=1:365
            if K_means(i)==1
               K_1(l,:) =PV_std(:,i);
               l=l+1;
            elseif K_means(i)==2
               K_2(m,:) =PV_std(:,i);
               m=m+1;
            elseif K_means(i)==3
               K_3(n,:) =PV_std(:,i);
               n=n+1;
            elseif K_means(i)==4
               K_4(o,:) =PV_std(:,i);
               o=o+1;
            end
        end
%         figure
% %         plot(K_1')
% %         hold on
%         plot(mean(K_1,1),'b--o')
%         hold on
% %         plot(K_2')
% %         hold on
%         plot(mean(K_2,1),'r--o')
%         hold on
% %         plot(K_3')
% %         hold on
%         plot(mean(K_3,1),'g--o');        
%         hold on
% %         plot(K_4')
% %         hold on
%         plot(mean(K_4,1),'y--o');  
        fourmean=[sum(mean(K_1,1)) sum(mean(K_2,1)) sum(mean(K_3,1)) sum(mean(K_4,1))];
        [Minmean4,INDEX4]= min(fourmean);
    D_start=0;
    for i=1:10
        for D=D_start+1:D_start+30
%            eval(['baseload',num2str(i),'(D,:)','=','NTL(1+(D-1)*96:96+(D-1)*96,i)',';'])%5-20
            eval(['ntlbase',num2str(i),'(D,:)','=','NTL(1+(D-1)*24:24+(D-1)*24,i)',';'])%29-80

        end
    end
    for i=1:10
%            eval(['baseload',num2str(i),'(D,:)','=','NTL(1+(D-1)*96:96+(D-1)*96,i)',';'])%5-20
            eval(['ntlbase',num2str(i),'(1:(D_start),:)','=','[]',';'])%29-80
    end    
    Max_1=min(ntlbase1);Max_2=min(ntlbase2);Max_3=min(ntlbase3);Max_4=min(ntlbase4);Max_5=min(ntlbase5);
    Max_6=min(ntlbase6);Max_7=min(ntlbase7);Max_8=min(ntlbase8);Max_9=min(ntlbase9);Max_10=min(ntlbase10);
    error1=sum(abs(ntlbase1-Max_1),2);error2=sum(abs(ntlbase2-Max_2),2);error3=sum(abs(ntlbase3-Max_3),2);error4=sum(abs(ntlbase4-Max_4),2);error5=sum(abs(ntlbase5-Max_5),2);
    error6=sum(abs(ntlbase6-Max_6),2);error7=sum(abs(ntlbase7-Max_7),2);error8=sum(abs(ntlbase8-Max_8),2);error9=sum(abs(ntlbase9-Max_9),2);error10=sum(abs(ntlbase10-Max_10),2);
    MAX_ALL=[Max_1;Max_2;Max_3;Max_4;Max_5;Max_6;Max_7;Max_8;Max_9;Max_10];
    error_all=[error1 error2 error3 error4 error5 error6 error7 error8 error9 error10];
    for i=1:10
    [Value_min Voting_min(i)]= min(error_all(:,i));
    end
    Day_mode=mode(Voting_min)+D_start;
%     Baseload_model=NTL((Day_mode-1)*24+1:24*Day_mode,:)';
%     [A_low,A_lowindex]=min(Baseload_model(:,1:6),[],2);

    A_low=MAX_ALL(:,3);
    %A_low=NTL(9+(Day_mode-1)*96:12+((Day_mode-1)*96),:)';
    repA=repmat(A_low,1,24);
    MAX_ALL_new=MAX_ALL-repA;
%     plot(MAX_ALL_new');
%    figure
    for i=1:10
        yy(i,:) = smooth(MAX_ALL_new(i,:),'lowess');
 %       plot(yy(i,:));
 %       hold on
    end
    %% no observation PV
%     prosumeryy=[yy(1:4,:);yy(6:7,:); yy(10,:)];
%     prosumeryy=mean(prosumeryy);
%     sumprosumeryy=abs((prosumeryy)/max(abs(prosumeryy)));
%     sumprosumeryy(1:32)=0;sumprosumeryy(71:96)=0;
%     PV_std(:,Day_mode)=sumprosumeryy';
    %%
    capacity_all=-yy./(PV_std(:,Day_mode)'+0.05);
    [max_capacity,index]=max(capacity_all(:,10:15),[],2);


    capacity=max_capacity;
    estimatepv=capacity*stdpv;%stdpv or 
%     figure
%     plot(estimatepv(3,1:D*24));
%     hold on
%     plot(onepv(3,1:D*24));
    estimatepv(10,:)=0;estimatepv(8,:)=0;estimatepv(5,:)=0;
    %onepv(10,:)=[];onepv(8,:)=[];onepv(5,:)=[];
    L_dis=abs(NTL+estimatepv');
    L_actual=abs(NTL+PV);
    L_dis(:,10)=L_actual(:,10);L_dis(:,8)=L_actual(:,8);L_dis(:,5)=L_actual(:,5);
    
    MAPE=(1/720)*sum(abs(L_dis(1:720,:)-L_actual(1:720,:))./(abs(L_actual(1:720,:))),1);
    MAPE_all=(1/720)*sum(abs(sum(L_dis(1:720,:),2)-sum(L_actual(1:720,:),2))./(abs(sum(L_actual(1:720,:),2))),1);
%     new_1=[];new_2=[];
%     for i=1:365
%         new_1=[new_1 estimatepv(:,36+(i-1)*96:66+(i-1)*96)];
%         new_2=[new_2 onepv(:,36+(i-1)*96:66+(i-1)*96)];
%     end
%     MAPE=(1/12775)*sum(abs(new_1-new_2)./(abs(new_2+0.1)),2);
    
    %MAPE_agg=1/35040*(sum(abs(sum(estimatepv,1)-sum(onepv,1))./abs(sum(onepv,1)+0.000001)));
    %MAPE_Day_mode_1=mean(sum(abs(estimatepv(:,(Day_mode-1)*96+1:Day_mode*96)-onepv(:,(Day_mode-1)*96+1:Day_mode*96)),2))./(sum(abs(onepv(:,(Day_mode-1)*96+1:Day_mode*96))+0.0001,2));    
    %MAPE_Day_mode_2=mean(sum(abs(estimatepv(:,(Day_mode-1)*96+1:Day_mode*96)-onepv(:,(Day_mode-1)*96+1:Day_mode*96)),2))./(sum(abs(onepv(:,(Day_mode-1)*96+1:Day_mode*96))+0.0001,2));
    
   % L_dis=abs(NTL+estimatepv');
%     Min_1=min(K_4);
%     avebaseload=mean(baseload,1);
%     a=1;
%     plot(oneweekdata,'DisplayName','oneweekdata')
%     hold on
%     z0=zeros(1,6720)';
%     plot(z0);
end