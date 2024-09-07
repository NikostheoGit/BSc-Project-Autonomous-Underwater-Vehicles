dt=0.8;
dl=0.8;
du=1.6;
d4=2;
l=1;
g=1;
met=1;
%diasthma [2,7]
metrhths3=0;
metrhths1=0;
metrhths2=0;
metrhths4=0;
f=0;
tsip=11;
metraw=0;
p88=1;
ari8mos=100;
Mesos=0;
Mesos1=0;
ari8mos1=200;

x(1)=1/2+2.5; %%Arxikoipoihsh polugwnou
y(1)=1/2+2.5;
x(2)=1+2.5;
y(2)=0.25+2.5;
x(3)=1.25+2.5;
y(3)=0.25+2.5;
x(4)=1.5+2.5;
y(4)=1/2+2.5;
x(5)=1.5+2.5;
y(5)=1.75+2.5;
x(6)=1.25+2.5;
y(6)=2+2.5;
x(7)=1+2.5;
y(7)=2+2.5;
x(8)=1/2+2.5;
y(8)=1.75+2.5;

for i=1:9
    if i==9
        c(i,1,1)=c(1,1);
        c(i,2,1)=c(1,2);
    else
        c(i,1,1)=x(i); %%Ftia3e pinaka me tis korufes tou polugwnou (Gia ti stigmh 1)
        c(i,2,1)=y(i);
    end
end


for m=2:50           %%Ftia3e 50 Plume
    for i=1:9
        
        c(i,1,m)=c(i,1,m-1) + (2*rand-1)*dt;
        c(i,2,m)=c(i,2,m-1) + (2*rand-1)*dt;
        
        if  c(i,1,m)> 7;
            c(i,1,m)=7;
        end
        if c(i,1,m)<2.25;
            c(i,1,m)=2.25;
        end
        
        if  c(i,2,m)> 7;
            c(i,2,m)=7;
        end
        if c(i,2,m)<2.25;
            c(i,2,m)=2.25;
        end
        
        if i==9
            c(i,1,m)=c(1,1,m);
            c(i,2,m)=c(1,2,m);
        end
    end
end

% 
%  figure(15)
% title(['Timestep 1 Initial Plume'])
%  axis([0 10 0 10]);
[a0,b0]=simple_polygon(8,c,1); %%Ftia3e to 1o polugwno
%  hold on
%  plot(a0,b0,'o-');
%  hold on

XY=[c(:,1,1),c(:,2,1)];
cent(1)=mean(c(:,1,1));
cent(2)=mean(c(:,2,1));
centgen(1,:)=[mean(c(:,1,1)),mean(c(:,2,1))];


plot(centgen(1,1),centgen(1,2),'x-');
hold on

Auv(1,:)=[centgen(1,1)-0.2,centgen(1,2)-0.1]; %%% Arxikopoihsh Auv Kentro
Auv1(1,:)=[c(5,1,1)+0.2,c(5,2,1)+0.2];   %%Arxikopoihsh Auv Korufwn
Auv1(2,:)=[c(8,1,1)-0.1,c(8,2,1)-0.1];
Auv1(3,:)=[c(1,1,1)+0.1,c(1,2,1)-0.1];

% rk = 2 + (7-2).*rand; %%3ekina apo to idio shmeio
% Auv(1,:)=[rk,rk]
% Auv1(1,:)=[rk,rk]
% Auv1(2,:)=[rk,rk]
% Auv1(3,:)=[rk,rk]
% 
% ra = 2 + (7-2).*rand; %%3ekina apo Random shmeio
% ya= 2 + (7-2).*rand
% Auv1(1,:)=[ra,ya]
% rb = 2 + (7-2).*rand; %%3ekina apo Random shmeio
% yb= 2 + (7-2).*rand
% Auv1(2,:)=[rb,yb]
% rc = 2 + (7-2).*rand; %%3ekina apo Random shmeio
% yc= 2 + (7-2).*rand
% Auv1(3,:)=[rc,yc]
% rd = 2 + (7-2).*rand; %%3ekina apo Random shmeio
% yd= 2 + (7-2).*rand
% Auv(1,:)=[rd,yd]


A1=Auv1(1,:)
A2=Auv1(2,:)
A3=Auv1(3,:)


[A,tmp1,tmp2]=EllipseDirectFit(XY); %%Bale thn Best Fit Ellipse sto arxiko polugwno(plume)
hold on
title(['Timestep 1 Initial Plume'])
per{1}=[tmp1,tmp2];


for i=10:10:50            %%Bres ta kentra barous stis stigmes ths anadyshs
    l=l+1;
    centgen(l,:)=[mean(c(:,1,i)),mean(c(:,2,i))]; %%Anadush
    p88=p88+1;
    
end

for i=5:10:50  %%Bres ta kentra barous sta mesodiasthmata
    f=f+1;
    centMes(f,:)=[mean(c(:,1,i)),mean(c(:,2,i))];%%Endiamesh Anadush
end

for g=10:10:50
    
    XY=[c(:,1,g),c(:,2,g)];
    cent=[mean(c(:,1,g)),mean(c(:,2,g))];
    
    
    plot(cent(1),cent(2),'x-');
    hold on
    
    [A,tmp1,tmp2]=EllipseDirectFit(XY);
    
    
    met=met+1; %ftia3e pinakes me ta shmeia twn ellipsewn.
    per{met}=[tmp1,tmp2];
end

for i=1:5
    [x,y]=simple_polygon(8,c,i*10);
    
     figure(ari8mos)

    for j=1:3 %%%% ALGORITHM 2 %%%%
        subplot(2,2,j)
        axis([0 10 0 10]);
        [lonc, latc]=circle(Auv1(j,1), Auv1(j,2), d4);
        hold on
        plot(Auv1(j,1), Auv1(j,2),'x-')
        hold on
        
        oldAuv1(j,:,i)=(Auv1(j,:));
        
        plot(x,y,'o-');
        hold on
        P2 = [x(:) y(:)];
        [F1 F2]=DataPolygon(P2);
        
        [loni, lati] = polyxpoly(lonc, latc,F1,F2);
        plot(loni,lati,'x')
        
        
        [rows columns] = size(loni);
        
        Sunt1(j,:,i)={loni ,lati};
    end
%         

[rows1 columns1]= size(Sunt1{1,1,i});
[rows2 columns2]= size(Sunt1{2,1,i});
[rows3 columns2]= size(Sunt1{3,1,i});



distance1=0
x1=0
y1=0
for i1=1:rows1
    for i2=1:rows2
        for i3=1:rows3
            
            
            distance_help=abs(Sunt1{2,1,i}(i2) -Sunt1{1,1,i}(i1))+abs(Sunt1{2,2,i}(i2)-Sunt1{1,2,i}(i1))+ abs(Sunt1{3,1,i}(i3) -Sunt1{1,1,i}(i1))+abs(Sunt1{3,2,i}(i3)-Sunt1{1,2,i}(i1))
            
            if distance_help>distance1
                
                distance1=distance_help
                x1=Sunt1{1,1,i}(i1)
                y1=Sunt1{1,2,i}(i1)
                
            end
            
            
        end
    end
end


distance2=0
x2=0
y2=0


for i2=1:rows2
        for i3=1:rows3
          
            
            distance_help2=abs(Sunt1{2,1,i}(i2) -x1)+abs(Sunt1{2,2,i}(i2)-y1)+ abs(Sunt1{3,1,i}(i3) -Sunt1{2,1,i}(i2))+abs(Sunt1{3,2,i}(i3)-Sunt1{2,2,i}(i2))
            
            if distance_help2>distance2
                
                distance2=distance_help2
                x2=Sunt1{2,1,i}(i2)
                y2=Sunt1{2,2,i}(i2)
               
		end 
       
        end
end



distance3=0
x3=0
y3=0



for i3=1:rows3

	distance_help3=abs(Sunt1{3,1,i}(i3) -x1)+abs(Sunt1{3,2,i}(i3)-y1)+ abs(Sunt1{3,1,i}(i3) -x2)+abs(Sunt1{3,2,i}(i3)-y2)


		if distance_help3>distance3
                
                distance3=distance_help3
                x3=Sunt1{3,1,i}(i3)
                y3=Sunt1{3,2,i}(i3)

       
        end
end

Auv1(1,:)=[x1 y1]
Auv1(2,:)=[x2 y2]
Auv1(3,:)=[x3 y3]

         for k=1:3
            if Auv1(k,1)==0 %CASE 3
                disp(['MPHKAN EDW'])
                
                [xunit,yunit]=circle( oldAuv1(k,1,i),oldAuv1(k,2,i),d4);
                for mol=1:9
                    for kol=1:9
                        R1(kol)=(c(kol,2,(i-1)*10+mol)-oldAuv1(k,2,i))/(c(kol,1,(i-1)*10+mol)-oldAuv1(k,1,i));
                        Mesos=Mesos+R1(kol);
                        
                    end
                end
                
                Mesos=Mesos/72;
                b1=oldAuv1(k,2)-Mesos*oldAuv1(k,1);
                [xout1,yout1] = linecirc(Mesos,b1,oldAuv1(k,1,i),oldAuv1(k,2,i),d4);
                
                p3=pdist2([xout1(1,1) yout1(1,1)],[centgen(i,1) centgen(i,2)]);
                p4=pdist2([xout1(1,2) yout1(1,2)],[centgen(i,1) centgen(i,2)]);
                
                if p3<p4
                    Auv1(k,:)=[abs(xout1(1,1)) abs(yout1(1,1))];
                end
                if p4>p3
                    Auv1(k,:)=[abs(xout1(1,2)) abs(yout1(1,2))];
                end
            end
         end

%             if rows>=2 %CASE 1
%                 Auv1(j-1,:)=[Sunt1{j-1,1,i}(m) Sunt1{j-1,2,i}(m)];
%                 Auv1(j,:)=[Sunt1{j,1,i}(n) Sunt1{j,2,i}(n)];
%             end
%             if rows==1 %CASE 2
%                 Auv1(j-1,:)=[Sunt1{j-1,1,i}(m) Sunt1{j-1,2,i}(m)];
%                 Auv1(j,:)=[Sunt1{j,1,i}(n) Sunt1{j,2,i}(n)];
%             end
%             title(['Auv ' num2str(j),'  TimeStep = ' num2str(i*10)])
%             if isempty(loni) %CASE 3
%                 [xunit,yunit]=circle(Auv1(j,1),Auv1(j,2),d4);
%                 for mol=1:9
%                     for kol=1:9
%                         R1(kol)=(c(kol,2,(i-1)*10+mol)-Auv1(j,2))/(c(kol,1,(i-1)*10+mol)-Auv1(j,1));
%                         Mesos=Mesos+R1(kol);
%                         
%                     end
%                 end
%                 
%                 Mesos=Mesos/72;
%                 b1=Auv1(j,2)-Mesos*Auv1(j,1);
%                 [xout1,yout1] = linecirc(Mesos,b1,Auv1(j,1),Auv1(j,2),d4);
%                 
%                 p3=pdist2([xout1(1,1) yout1(1,1)],[centgen(i,1) centgen(i,2)]);
%                 p4=pdist2([xout1(1,2) yout1(1,2)],[centgen(i,1) centgen(i,2)]);
%                 
%                 if p3<p4
%                     Auv1(j,:)=[xout1(1,1) yout1(1,1)];
%                 end
%                 if p4>p3
%                     Auv1(j,:)=[xout1(1,2) yout1(1,2)];
%                 end
%                 
%                 

    
    %%%%%%%%%---Telos Algori8mos 2---%%%%
   
    
    ari8mos=ari8mos+1;
    
    %%%%%%%---ALGORITHMOS 1---%%%%%
    F(i)=pdist2(Auv(i,:),centgen(i+1,:));
    
    %%%%%%%%----------Case 1------------%%%%%%%%%%%
    if dl<F(i) && F(i)<du
        figure(i)
        fig1=scatter(Auv(i,1),Auv(i,2),'g','filled');
        hold on
        plot(x,y,'o-');
        hold on
        XY=[c(:,1,i*10),c(:,2,i*10)];
        [A,tmp1,tmp2]=EllipseDirectFit(XY);
        hold on
        metraw=1;
        metrhths1=metrhths1+1;
        Auv(i+1,:)=centgen(i+1,:);
        
        figure(i)
        fig2=scatter(Auv(i+1,1),Auv(i+1,2),'r','filled');
        title(['Case 1 Timestep = ' num2str(i*4)])
        legend([fig1 fig2],'Arxikh 8esh Auv','Telikh 8esh Auv')
        hold on
        
    end
    %%%%%%%%----------Case 2------------%%%%%%%%%%%
    if F(i)<dl
        
        matrixper= per{1,i+1};
        for j=1:size(matrixper)
            G(j)=pdist2(Auv(i,:),matrixper(j,:));
        end
        minimum1=min(G);
        [skoros1] = min(G(:)) ;
        [bres1] = ind2sub(size(G),find(G==skoros1));
        
        E=i*10;
        for k=1:9
            Q(k)=pdist2(Auv(i,:),c(k,:,E));
        end
        minimum2=min(Q);
        [skoros2] = min(Q(:)) ;
        [bres2] = ind2sub(size(G),find(Q==skoros2));
        
        minmumALL=min(minimum1,minimum2);
        
        if dl<=minmumALL && minmumALL<=du
            figure(i)
            
            fig3=scatter(Auv(i,1),Auv(i,2),'g','filled');
            hold on
            XY=[c(:,1,i*10),c(:,2,i*10)];
            [A,tmp1,tmp2]=EllipseDirectFit(XY);
            hold on
            plot(x,y,'o-');
            hold on
            metraw=1;
            metrhths2=metrhths2+1;
            if minimum1<minimum2
                Auv(i+1,:)=matrixper(bres1,:);
                figure(i)
                fig4=scatter(Auv(i+1,1),Auv(i+1,2),'y','filled');
                hold on
                Auv(i+1,:)=centgen(i+1,:);
                figure(i)
                fig5=scatter(Auv(i+1,1),Auv(i+1,2),'r','filled');
                hold on
                title(['Case 2 (Ellipsh) Timestep = ' num2str(i*4)])
                legend([fig3 fig4 fig5],'Arxikh 8esh Auv','Endiamesh 8esh Auv','Telikh 8esh Auv')
            else
               
                Auv(i+1,:)=c(bres2,:,E);
                figure(i)
                fig6=scatter(Auv(i+1,1),Auv(i+1,2),'y','filled');
                hold on
                Auv(i+1,:)=centgen(i+1,:);
                figure(i)
                fig7=scatter(Auv(i+1,1),Auv(i+1,2),'r','filled');
                hold on
                title(['Case 2 (Korufes) Timestep = ' num2str(i*4)])
                legend([fig3 fig6 fig7],'Arxikh 8esh Auv','Endiamesh 8esh Auv','Telikh 8esh Auv')
            end
        end
        
        
    end
    
    %%%%%%%%----------Case 3------------%%%%%%%%%%%
    if F(i)<dl
        if metraw==0
            
            figure(i)
            axis([0 10 0 10]);
            fig8=scatter(Auv(i,1),Auv(i,2),'g','filled');
            hold on
            [a,k0,l0]=Superellipse(centgen(i+1,:),Auv(i,:));
            Auv(i+1,:)=[k0(7) l0(7)];
            fig9=scatter(Auv(i+1,1),Auv(i+1,2),'y','filled');
            hold on
            XY=[c(:,1,i*10),c(:,2,i*10)];
            [A,tmp1,tmp2]=EllipseDirectFit(XY);
            hold on
            
            plot(x,y,'o-');
            hold on
            Auv(i+1,:)=centgen(i+1,:);
            fig10=scatter(Auv(i+1,1),Auv(i+1,2),'r','filled');
            hold on
            title(['Case 3 Timestep = ' num2str(i*4)])
            legend([fig8 fig9 fig10],'Arxikh 8esh Auv','Endiamesh 8esh Auv','Telikh 8esh Auv')
            metrhths3=metrhths3+1;
            metraw=1;
        end
    end
    
    %%%%%%%%----------Case 4------------%%%%%%%%%%%
    if  F(i)> du
        metrhths4=metrhths4+1;
        [xunit,yunit]=circle(Auv(i,1),Auv(i,2),d4);
        figure(i)
        
        
        fig11=scatter(Auv(i,1),Auv(i,2),'g','filled');
        hold on
        plot(x,y,'o-');
        hold on
        XY=[c(:,1,i*10),c(:,2,i*10)];
        [A,tmp1,tmp2]=EllipseDirectFit(XY);
        hold on
        R=(centMes(i,2)-Auv(i,2))/(centMes(i,1)-Auv(i,1));
        
        b=Auv(i,2)-R*Auv(i,1);
        
        [xout,yout] = linecirc(R,b,Auv(i,1),Auv(i,2),d4);
        
        p1=pdist2([xout(1,1) yout(1,1)],[centMes(i,1) centMes(i,2)]);
        p2=pdist2([xout(1,2) yout(1,2)],[centMes(i,1) centMes(i,2)]);
        
        if p1<p2
            Auv(i+1,:)=[xout(1,1) yout(1,1)];
             figure(i)
            fig12=scatter(Auv(i+1,1),Auv(i+1,2),'r','filled');
            hold on
            title(['Case 4 Timestep = ' num2str(i*4)])
            legend([fig11 fig12],'Arxikh 8esh Auv','Telikh 8esh Auv')
        end
        if p1>p2
            Auv(i+1,:)=[xout(1,2) yout(1,2)];
            figure(i)
            fig13=scatter(Auv(i+1,1),Auv(i+1,2),'r','filled');
            hold on
            title(['Case 4 Timestep = ' num2str(i*10)])
            legend([fig11 fig13],'Arxikh 8esh Auv','Telikh 8esh Auv')
        end
    end

    
    metraw=0;
    figure(6)
      axis([0 10 0 10]);
     plot(a0,b0,'o-');
     hold on
   
      fig3=scatter(Auv(1,2),Auv(1,2),'b','filled');
      hold on
      fig0=scatter(A1(1),A1(2),'g','filled');
      hold on
      fig1=scatter(A2(1),A2(2),'g','filled');
      hold on
      fig2=scatter(A3(1),A3(2),'g','filled');
      title(['Timestep 1 Initial Plume'])
      legend([fig1 fig3],'Arxiko perifereia OldAuv1','Arxiko Auv kentro')    
      saveas(figure(6),sprintf('arxiko.png'))
    h=figure(ari8mos1)
    plot(x,y,'o-');
   
    hold on
    axis([0 10 0 10]);
   % fig14= scatter(oldAuv1(:,1,i),oldAuv1(:,2,i),'g','filled');
%    labels = {'OldAuv1' 'OldAuv2' 'OldAuv3'};
%    labelpoints(oldAuv1(:,1,i),oldAuv1(:,2,i), labels);
%    hold on
%    labels1 = {'NewAuv1' 'NewAuv2' 'NewAuv3'};
%    labelpoints(Auv1(:,1,i),Auv1(:,2,i), labels1);
    hold on
    fig15=scatter(Auv1(:,1),Auv1(:,2),'r','filled');
    hold on
%     fig16=scatter(Auv(i,1),Auv(i,2),'b','filled');
%     hold on
    fig17=scatter(Auv(i+1,1),Auv(i+1,2),'b','filled');
    hold on
    title(['All Timestep = ' num2str(i*4)])
    l1= legend([fig15 fig17],'Teliko perifereia Auv1','Teliko Auv kentro')
    ari8mos1=ari8mos1+1;
%     saveas(h,sprintf('FIG%d.png',ari8mos1)); % will create FIG1, FIG2,...
%     saveas(h,sprintf('FIG%d.fig',ari8mos1)); % will create FIG1, FIG2,...
    
end
close(100)
close(101)
close(102)
close(103)
close(104)
close(201)
close(202)
close(203)
close(204)







