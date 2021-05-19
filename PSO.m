function [sbestpos,sbestval,objit] = PSO(lb_lamdas,ub_lamdas,d_lamdas,lb_nonIt3,ub_nonIt3,d_nonIt3,lb_It4,ub_It4,d_It4,ub_gamma,lb_gamma,d_gamma,ssize,w,c1,c2,iter_number)
    %PSO Summary of this function goes here
    %   Detailed explanation goes here

    
    
    d=d_lamdas+d_nonIt3+d_It4+d_gamma; % number of variables
    

    %%%%% first values of swarm members are selected as heuristically %%%%%
    swarm=[unifrnd(lb_lamdas,ub_lamdas,[ssize,d_lamdas]),...
        unifrnd(lb_nonIt3,ub_nonIt3,[ssize,d_nonIt3]),...
        unifrnd(lb_It4,ub_It4,[ssize,d_It4]),...
        unifrnd(lb_gamma,ub_gamma,[ssize,d_gamma])];   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    obj=zeros(ssize,1);          % array defined for objective values
   
    % max boundary of velocity of change of value of particles
    % max boundary value of velocity selected as half of range by 
    % literature, generally.
    
    vmax_lamdas=(ub_lamdas-lb_lamdas)/64;      
    vmax_nonIt3=(ub_nonIt3-lb_nonIt3)/64;
    vmax_It4=(ub_It4-lb_It4)/64;
    vmax_gamma=(ub_gamma-lb_gamma)/64;
   
    % first objective values are calculated for all swarm. 
    for i=1:ssize
       obj(i)=DynamicProcess(swarm(i,:)); 

    end

    % first velocity value can be selected zero or very little value
    % generally.(literature information)
    
    velocity=zeros(ssize,d)+0.001.*ones(ssize,d);
    
    
    pbestpos=swarm; % first best position is equal to swarm
    pbestval=obj;   % first best objective value equal to objective value 

    
    [sbestval idx]=min(obj); 
    sbestpos=swarm(idx,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%% ITERATION of PSO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for iter=1:iter_number
        
        % velocity determine as heuristically
        for i=1:ssize
            velocity(i,:)=w.*velocity(i,:)+c1.*unifrnd(0,1).*...
                (pbestpos(i,:)-swarm(i,:))+c2.*unifrnd(0,1).*...
                (sbestpos-swarm(i,:));

        end

        % boundary validation of determined velocity values  
        for i=1:ssize
            
            for j=1:d_lamdas
                if velocity(i,j)>vmax_lamdas
                    velocity(i,j)=vmax_lamdas;
                elseif velocity(i,j)<-vmax_lamdas
                    velocity(i,j)=-vmax_lamdas;     
                end

            end
            
            
            if velocity(i,4)>vmax_nonIt3
                    velocity(i,4)=vmax_nonIt3;
            elseif velocity(i,4)<-vmax_nonIt3
                    velocity(i,4)=-vmax_nonIt3;     
            end
            
            
            if velocity(i,5)>vmax_It4
                    velocity(i,5)=vmax_It4;
            elseif velocity(i,5)<-vmax_It4
                    velocity(i,5)=-vmax_It4;     
            end
            
            if velocity(i,6)>vmax_gamma
                    velocity(i,6)=vmax_gamma;
            elseif velocity(i,6)<-vmax_gamma
                    velocity(i,6)=-vmax_gamma;     
            end


        end
        
        swarm=swarm+velocity; % update of swarm
        
        % boundary validation of updated swarm 
        for i=1:ssize
            for j=1:d_lamdas
                if swarm(i,j)>ub_lamdas
                    swarm(i,j)=ub_lamdas;
                elseif swarm(i,j)<lb_lamdas
                    swarm(i,j)=lb_lamdas;     
                end
                
            end
            if swarm(i,4)>ub_nonIt3
                    swarm(i,4)=ub_nonIt3;
            elseif swarm(i,4)<lb_nonIt3
                    swarm(i,4)=lb_nonIt3;     
            end
            if swarm(i,5)>ub_It4
                    swarm(i,5)=ub_It4;
            elseif swarm(i,5)<lb_It4
                    swarm(i,5)=lb_It4;     
            end
            
            if swarm(i,6)>ub_gamma
                    swarm(i,6)=ub_gamma;
            elseif swarm(i,6)<lb_gamma
                    swarm(i,6)=lb_gamma;     
            end
            
        end
        
        % 'iter'th objective value calculation 
        for i=1:ssize
           obj(i)=DynamicProcess(swarm(i,:)); 

        end

        %best partical update
        for i=1:ssize
            if obj(i)<pbestval(i)
                pbestval(i)=obj(i);
                pbestpos(i,:)=swarm(i,:);
            end
        end
            
        % best swarm update
        if min(obj)<sbestval
            [sbestval idx]=min(obj);
             sbestpos=swarm(idx,:);
        end
        objit(iter)=sbestval;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end