%--- Start ---%

    % Tidy up and make the window clean
        clc;
        clear;

%--- Parameter Sections ---%
    pitchY=0.82;
    width_total=3.63;
    width_walk=0.64;
    widthX=(width_total-width_walk)/2;
    front_space=3;
    T_i = 0.1;
    pitch_total = 24*pitchY;
    U = 10;
    R=0.2;
    smallx = 0.1;
    PS= 0.5;
    A = 0.001;
    B = 0.08;
    K= 2.4;
    k = 1.2;
    count = 0;

    
    % Main parameters
            crowdsize = 144; 
            timestep = 1;
            total_t = 2000;
            timelimit =1:timestep:total_t;
            luggage_percentage = 80;
            luggage_time = 7.2;
            narrow_coef = 0.2;
            
            %Sex% distribution
            men = 65;
            woman = 35;
            % Age% distribution accumulated. Change only the number.
            age20= 17;
            age30= 23.5 + age20;
            age40= 23.5 + age30;
            age50= 19 + age40;
            age60= 20 + age50;
            age70= 4 + age60;
            
           % Mean weight for men (kg)
           mean20m = 80.7;
            sd20m = 9.4;
           mean30m = 80.4;
            sd30m = 13.1;
           mean40m = 86.2;
            sd40m = 13.5;
           mean50m = 87.2;
            sd50m = 14.3;
           mean60m = 81.1;
            sd60m = 9.8;
           mean70m = 77.3;
            sd70m = 9.1;
            
            % Speed for men (cm/s)
            meanv20m = 139.3;
            sdv20m = 15.3;
           meanv30m = 145.8;
            sdv30m = 9.4;
           meanv40m = 146.2;
            sdv40m = 16.4;
           meanv50m = 139.3;
            sdv50m = 22.9;
           meanv60m = 135.9;
            sdv60m = 20.5;
           meanv70m = 133;
            sdv70m = 19.6;
            
            % Mean weight for women (kg)
           mean20w = 59;
            sd20w = 6.5;
           mean30w = 66.4;
            sd30w = 18.9;
           mean40w = 63.2;
            sd40w = 13.2;
           mean50w = 64.5;
            sd50w = 11.4;
           mean60w = 63.4;
            sd60w = 9.7;
           mean70w = 59.3;
            sd70w = 8.6;

            % Speed for women (cm/s)
            meanv20w = 140.7;
            sdv20w = 17.5;
           meanv30w = 141.5;
            sdv30w = 12.7;
           meanv40w = 139.1;
            sdv40w = 15.8;
           meanv50w = 139.5;
            sdv50w = 15.1;
           meanv60w = 129.6;
            sdv60w = 21.3;
           meanv70w = 127.2;
            sdv70w = 21.1;
            
    % Control panel (0 = off and 1=on)
    % status has to have only 1 turned on. e.g. can't turn on Wilma and RP
    % on at the same time.
    
            processT = 9;
            seat1_interference = 0;
            seat2_interference = 0;
            row_max = ceil(crowdsize/6);
    
	% All about Wilma Method
            wilma_status = 0;
            limitcallW = 300;
            wilma_forcalibrate = 0;     %Should be 0 for normal cases
            
            num_groupW1 = row_max*2;
            num_groupW2 = row_max*2;
            num_groupW3 = row_max*2;
            
            start_groupW1 = 0;

                if num_groupW1*processT > limitcallW && wilma_forcalibrate == 0
                    start_groupW2 = start_groupW1 + (num_groupW1*processT);
                else
                    start_groupW2 = start_groupW1 + limitcallW;
                end

                if num_groupW2*processT > limitcallW && wilma_forcalibrate == 0
                    start_groupW3 = start_groupW2 + (num_groupW2*processT);
                else
                    start_groupW3 = start_groupW2 + limitcallW;  
                end
                
                if wilma_forcalibrate == 1
                     start_groupW1 = 0;
                     start_groupW2 = start_groupW1 + (2*row_max*processT) + processT;
                     start_groupW3 = start_groupW2 + (2*row_max*processT) + processT;
                end

            GW12_difT = start_groupW2-start_groupW1;
            GW23_difT = start_groupW3-start_groupW2;
            
    % All about Block method
            block_status = 0;
            num_groupB1 = crowdsize/3;
            num_groupB2 = crowdsize/3;
            num_groupB3 = crowdsize/3;
            limitcallB = 300;
            
            start_groupB1 = 0;
            start_groupB2 = start_groupB1 + (num_groupB1*processT);
            start_groupB3 = start_groupB2 + (num_groupB2*processT);

	% All about Random
            random_status = 1;
            num_groupR1 = crowdsize/6;
            num_groupR2 = crowdsize/6;
            num_groupR3 = crowdsize/6;
            num_groupR4 = crowdsize/6;
            num_groupR5 = crowdsize/6;
            num_groupR6 = crowdsize/6;
            limitcallR = 300;
            
             start_groupR1 = 0;
             start_groupR2 = start_groupR1 + (num_groupR1*processT);
             start_groupR3 = start_groupR2 + (num_groupR2*processT);
             start_groupR4 = start_groupR3 + (num_groupR3*processT);
             start_groupR5 = start_groupR4 + (num_groupR4*processT);
             start_groupR6 = start_groupR5 + (num_groupR5*processT);
            
    % All about Reversed pyramid
        
            RP_status = 0;      %Turn on and off 0/1

            group3plus_percent = 30;
            group2_percent = 30;
            row3plus_require = ceil(group3plus_percent/100*crowdsize/6);
            row2_require = ceil(group2_percent/100*crowdsize/4);
            
           
            % This is an unfinished algo to find the optimal row arrangement. It works only to some extent with lot of uncovered cases. 
                    controlrow23 = row_max - row3plus_require+1;
                    controlrow12 = controlrow23 - row2_require;
                    controlrow45 = controlrow23 + ceil(0.25*(row_max-controlrow23));
                    controlrow34 = controlrow45-1;

                    if controlrow12 > controlrow23 || controlrow12 == controlrow23
                        controlrow23 = controlrow12 + 1;
                    end
                    if controlrow23 > controlrow34 || controlrow23 == controlrow34
                        controlrow34 = controlrow23 + 1;
                    end
                    if controlrow34 > controlrow45 || controlrow34 == controlrow45
                    controlrow45 = controlrow34 + 1;
                    end
            
            
            num_group1 = (row_max-controlrow45+1)*6;
            num_group3 = (controlrow23-1)*2 + (controlrow12-1)*2;
            num_group2 = crowdsize - num_group1 - num_group3;
             
            start_group1 = 0;
            limitcall = 300;
            
                if  num_group1*processT > limitcall
                    start_group2 = num_group1*processT;
                else
                     start_group2 = start_group1 + limitcall;
                end

                if  num_group2*processT >  limitcall
                    start_group3 = start_group2   +   (num_group2*processT);
                else
                    start_group3 = start_group2 + limitcall;
                end
            
            G12_difT = start_group2-start_group1;
            G23_difT = start_group3-start_group2;



%--- Managing Arrays and pre-allocation Section for optimal speed---%

    Coordinates = zeros(2,1,crowdsize,length(timelimit));
    Velocities = zeros(2,1,crowdsize,length(timelimit));
    Accelerations = zeros(2,1,crowdsize,length(timelimit));
    Forces = zeros(2,1,crowdsize,length(timelimit));
    M = zeros(crowdsize,2);
    mass = zeros(crowdsize,1);
    mass_agent = zeros(crowdsize,1);
    Desireddirections = zeros(2,1,crowdsize,length(timelimit));
    Desiredspeed = zeros(crowdsize,length(timelimit));
    speed_agent = zeros(crowdsize,1);
    Maxspeed = zeros(crowdsize,length(timelimit));
    X = zeros(crowdsize,length(timelimit));
    Y = zeros(crowdsize,length(timelimit));
    diameter = linspace(0.5,0.7,crowdsize);
    inplane = zeros(crowdsize,1);
    inseat = zeros(crowdsize,1);
    blocknmoving = zeros(crowdsize,1);
    release_time =  zeros(crowdsize,1);
    luggage =  zeros(crowdsize,1);
    start_lug = zeros(crowdsize,1);
    stop_lug = zeros(crowdsize,1);
    waitlug = zeros(crowdsize,1);
    timetoseat = zeros(crowdsize,1);

    
    %Validations pre-allocation
    fail = zeros(crowdsize,1);
    fail2 = zeros(crowdsize,1);
    
%--- Randomize the seat ---%
seat =  randperm(crowdsize,crowdsize);

%--- Life is a random thing ---%
life = rand(crowdsize,1);
life2 = rand(crowdsize,1);


% Determining the row of the seat
for i=1:crowdsize
    row = ceil(seat/6);
    rowy = (24-row)*pitchY+0.5*pitchY;
end
    
% Managing the luggage (yes or no). Randomize it and then determine
 lug= randperm(ceil(crowdsize/(luggage_percentage/100)),crowdsize);
for i=1:crowdsize
  if lug(i) > crowdsize
      luggage(i) = 0;
  else 
      luggage(i) = 1;
  end
end

    
%--- Map Section ---%

    % Coordinate of the x-axis
    planx = [widthX,0,0,  widthX,0,0,  widthX,0,0,  widthX,0,0,  widthX,0,0,  widthX,0,0,  widthX,0,0,  widthX,0,0,  widthX,0,0,  widthX,0,0....
        ,widthX,0,0,  widthX,0,0,  widthX,0,0,  widthX,0,0,  widthX,0,0,  widthX,0,0,  widthX,0,0,  widthX,0,0,  widthX,0,0,  widthX,0,0 ....
        ,widthX,0,0,  widthX,0,0,  widthX,0,0,  widthX,0,0,  widthX,0 ];
    % Coordinate of the y-axis
    plany = [0,0,pitchY*1, pitchY*1,pitchY*1,pitchY*2, pitchY*2,pitchY*2,pitchY*3, pitchY*3,pitchY*3,pitchY*4, pitchY*4,pitchY*4,pitchY*5....
        pitchY*5,pitchY*5,pitchY*6,pitchY*6,pitchY*6,pitchY*7,pitchY*7,pitchY*7,pitchY*8,pitchY*8,pitchY*8,pitchY*9,pitchY*9,pitchY*9,pitchY*10....
        ,pitchY*10,pitchY*10,pitchY*11,pitchY*11,pitchY*11,pitchY*12,pitchY*12,pitchY*12,pitchY*13,pitchY*13,pitchY*13,pitchY*14,pitchY*14,pitchY*14,pitchY*15....
        ,pitchY*15,pitchY*15,pitchY*16,pitchY*16,pitchY*16,pitchY*17,pitchY*17,pitchY*17,pitchY*18,pitchY*18,pitchY*18,pitchY*19,pitchY*19,pitchY*19,pitchY*20....
        ,pitchY*20,pitchY*20,pitchY*21,pitchY*21,pitchY*21,pitchY*22,pitchY*22,pitchY*22,pitchY*23,pitchY*23,pitchY*23,pitchY*24,pitchY*24,pitchY*24];
    % Trick to plot the right-hand seats
	planrx = width_total - planx;
    % Top section of the walkway
	planTopx = [0,widthX+width_walk,widthX+width_walk];
        planTopy = [pitch_total+front_space,pitch_total+front_space,pitch_total];
    % Bottom section of the walkway to fill the border
	planBtmx = [widthX,widthX+width_walk];
        planBtmy = [0,0];
    % Customer spawning section
    spawnx = -5;
        planSpawnx = [0 ,spawnx ,spawnx, widthX, widthX];
        planSpawny = [pitch_total+front_space,pitch_total+front_space,pitch_total+1,pitch_total+1,pitch_total];

        randomSpawnx = randi([0 100],crowdsize,1)/100;
        randomSpawny = randi([-100 100],crowdsize,1)/100;


%--- Agent info section ---%
        
sex = randperm(crowdsize,crowdsize);
age = randperm(crowdsize,crowdsize);

% 1 for men and 0 for woman
% 2 for age20s, 3 for 30s and so on
for i=1:crowdsize
  % Determine sex
  if sex(i) > men*crowdsize/100
      sex(i) = 0;
  else 
      sex(i) = 1;
  end
  
  % Determine age
  if age(i) < age20*crowdsize/100
      age(i) = 2;
  elseif age(i) < age30*crowdsize/100 && age(i) >= age20*crowdsize/100
      age(i) = 2;
  elseif age(i) < age40 *crowdsize/100 && age(i) >= age30*crowdsize/100
      age(i) = 3;
  elseif age(i) < age50*crowdsize/100  && age(i) >= age40*crowdsize/100
      age(i) = 4;
  elseif age(i) < age60*crowdsize/100 && age(i) >= age50*crowdsize/100
      age(i) = 5;    
  elseif age(i) < age70*crowdsize/100 && age(i) >= age60*crowdsize/100
      age(i) = 6;
  elseif age(i) >= age70*crowdsize/100   
      age(i) = 7;  
  end
  
    if sex(i) == 1 && age(i) == 2
        mass_agent(i) = normrnd(mean20m,sd20m);
        speed_agent(i) = normrnd(meanv20m,sdv20m)/100 * narrow_coef;
    elseif sex(i) == 1 && age(i) == 3
        mass_agent(i) = normrnd(mean30m,sd30m);
        speed_agent(i) = normrnd(meanv30m,sdv30m)/100 * narrow_coef;
    elseif sex(i) == 1 && age(i) == 4
        mass_agent(i) = normrnd(mean40m,sd40m);
        speed_agent(i) = normrnd(meanv40m,sdv40m)/100* narrow_coef;
    elseif sex(i) == 1 && age(i) == 5
        mass_agent(i) = normrnd(mean50m,sd50m);
        speed_agent(i) = normrnd(meanv50m,sdv50m)/100* narrow_coef;
    elseif sex(i) == 1 && age(i) == 6
        mass_agent(i) = normrnd(mean60m,sd60m);
        speed_agent(i) = normrnd(meanv60m,sdv60m)/100* narrow_coef;
    elseif sex(i) == 1 && age(i) == 7
        mass_agent(i) = normrnd(mean70m,sd70m);
        speed_agent(i) = normrnd(meanv70m,sdv70m)/100* narrow_coef;
    end
    
    if sex(i) == 0 && age(i) == 2
        mass_agent(i) = normrnd(mean20w,sd20w);
        speed_agent(i) = normrnd(meanv20w,sdv20w)/100* narrow_coef;
    elseif sex(i) == 0 && age(i) == 3
        mass_agent(i) = normrnd(mean30w,sd30w);
        speed_agent(i) = normrnd(meanv30w,sdv30w)/100* narrow_coef;
    elseif sex(i) == 0 && age(i) == 4
        mass_agent(i) = normrnd(mean40w,sd40w);
        speed_agent(i) = normrnd(meanv40w,sdv40w)/100* narrow_coef;
    elseif sex(i) == 0 && age(i) == 5
        mass_agent(i) = normrnd(mean50w,sd50w);
        speed_agent(i) = normrnd(meanv50w,sdv50w)/100* narrow_coef;
    elseif sex(i) == 0 && age(i) == 6
        mass_agent(i) = normrnd(mean60w,sd60w);
        speed_agent(i) = normrnd(meanv60w,sdv60w)/100* narrow_coef;
    elseif sex(i) == 0 && age(i) == 7
        mass_agent(i) = normrnd(mean70w,sd70w);
        speed_agent(i) = normrnd(meanv70w,sdv70w)/100* narrow_coef;
    end
end
        
      
%--- Initialization Section (t=1) ---%
for i = 1:crowdsize
    mass(i,1) = mass_agent(i);
    Velocities(:,:,i,1) = [0;0];
    Coordinates(:,:,i,1) = [spawnx+2.5-randomSpawnx(i);pitch_total+1+0.5*(front_space-1)+(0.1*randomSpawny(i)*(front_space-1))];
    %Coordinates(:,:,i,1) = [-4;pitch_total+2];
    Desiredspeed(i,1) = speed_agent(i);
    Maxspeed(i,1) = Desiredspeed(i,1)*1.3;
    Desireddirections(:,:,i,1) = [0;0];
    if norm(Velocities(:,:,i,1)) > Maxspeed(i,1) % limit velocity to max speed
        arg = Velocities(2,1,i,1)/Velocities(1,1,i,1);
        if Velocities(1,1,i,t) < 0
            Velocities(:,:,i,t) = [-Maxspeed(i,t)*cos(arg);Maxspeed(i,t)*sin(arg)];
        else
            Velocities(:,:,i,t) = [Maxspeed(i,t)*cos(arg);Maxspeed(i,t)*sin(arg)];
        end
    end
    Forces(:,:,i,1) = Forces(:,:,i,1) + (1/T_i).*(Desiredspeed(i,1).*Desireddirections(:,:,i,1) - Velocities(:,:,i,1)); % Acceleration/deceleration force
    Accelerations(:,:,i,1) = Forces(:,:,i,1)./mass(i,1);
end

%Queing system
for i = 1:crowdsize
                AAgent = find(seat==(row(i)*6-5));
                BAgent = find(seat==(row(i)*6-4));
                CAgent = find(seat==(row(i)*6-3));
                DAgent = find(seat==(row(i)*6-2));
                EAgent = find(seat==(row(i)*6-1));
                FAgent = find(seat==(row(i)*6-0));
                
         if block_status ==1
             if row(i) > 2/3*row_max
                 release_time(i) = start_groupB1 + (randperm(num_groupB1,1) * processT) +ceil(limitcallB*rand);
             elseif row(i) < 1/3*row_max
                 release_time(i) = start_groupB2 + (randperm(num_groupB2,1) * processT) +ceil(limitcallB*rand);
             elseif row(i) >= 1/3*row_max && row(i) <= 2/3*row_max
                 release_time(i) = start_groupB3 + (randperm(num_groupB3,1) * processT) +ceil(limitcallB*rand);
             end         
         end
         
        if random_status ==1

            if mod(i,6) == 0
            release_time(i) = start_groupR1 + (randperm(num_groupR1,1) * processT) +ceil(limitcallR*rand);
            end
            if mod(i,6) == 1
            release_time(i) = start_groupR2 + (randperm(num_groupR2,1) * processT) +ceil(limitcallR*rand);
            end
            if mod(i,6) == 2
            release_time(i) = start_groupR3 + (randperm(num_groupR3,1) * processT) +ceil(limitcallR*rand);
            end
            if mod(i,6) == 3
            release_time(i) = start_groupR4 + (randperm(num_groupR4,1) * processT) +ceil(limitcallR*rand);
            end
            if mod(i,6) == 4
            release_time(i) = start_groupR5 + (randperm(num_groupR5,1) * processT) +ceil(limitcallR*rand);
            end
            if mod(i,6) == 5
            release_time(i) = start_groupR6 + (randperm(num_groupR6,1) * processT) +ceil(limitcallR*rand);
            end
        end
        
        if and(wilma_status == 1,wilma_forcalibrate == 0)
                    release_time(AAgent) = start_groupW1 + ceil(processT * randperm(num_groupW1,1)+ceil(limitcallW*rand));
                    release_time(FAgent) = start_groupW1 + ceil(processT * randperm(num_groupW1,1)+ceil(limitcallW*rand));
                    release_time(BAgent) = start_groupW2 + ceil(processT * randperm(num_groupW2,1)+ceil(limitcallW*rand));
                    release_time(EAgent) = start_groupW2 + ceil(processT * randperm(num_groupW2,1)+ceil(limitcallW*rand));
                    
                    % The last release time for last group has to be
                    % controlled, so that one agent is always released at
                    % the latest control time.
                    
                     if row(i) == ceil(row_max/4)
                        release_time(CAgent) = start_groupW3 + ceil(processT * (num_groupW3-2)) +ceil(limitcallW*rand); 
                     else
                        release_time(CAgent) = start_groupW3 + ceil(processT * randperm(num_groupW3,1)+ceil(limitcallW*rand));
                     end
                    
                     if row(i) == ceil(3*row_max/4)
                        release_time(DAgent) = start_groupW3 + ceil(processT * (num_groupW3-1))+ceil(limitcallW*rand);
                     else
                        release_time(DAgent) = start_groupW3 + ceil(processT * randperm(num_groupW3,1)+ceil(limitcallW*rand));
                     end
        end
        
        if and(wilma_status == 1  , wilma_forcalibrate == 1)
             release_time(AAgent) = start_groupW1 + (row_max - row(i))*processT;
             release_time(FAgent) = start_groupW1 + (row_max * processT) + ((row_max - row(i))*processT);
             release_time(BAgent) = start_groupW2 + (row_max - row(i))*processT;
             release_time(EAgent) = start_groupW2 + (row_max * processT) + ((row_max - row(i))*processT);
             release_time(CAgent) = start_groupW3 + (row_max - row(i))*processT;
             release_time(DAgent) = start_groupW3 + (row_max * processT) + ((row_max - row(i))*processT);
        end
        
        
        
        if RP_status == 1
                    if row(i) >= controlrow45
                        AAgent45 = find(seat==(row(i)*6-5));
                        BAgent45 = find(seat==(row(i)*6-4));
                        CAgent45 = find(seat==(row(i)*6-3));
                        DAgent45 = find(seat==(row(i)*6-2));
                        EAgent45 = find(seat==(row(i)*6-1));
                        FAgent45 = find(seat==(row(i)*6-0));

                        release_time(AAgent45) = start_group1 + (processT * randperm(num_group1,1)+ceil(limitcall*rand));
                        release_time(FAgent45) = start_group1 + (processT * randperm(num_group1,1)+ceil(limitcall*rand));
                        release_time(BAgent45) = start_group1 + (processT * randperm(num_group1,1)+ceil(limitcall*rand));
                        release_time(EAgent45) = start_group1 + (processT * randperm(num_group1,1)+ceil(limitcall*rand));
                        release_time(CAgent45) = start_group1 + (processT * randperm(num_group1,1)+ceil(limitcall*rand));
                        release_time(DAgent45) = start_group1 + (processT * randperm(num_group1,1)+ceil(limitcall*rand));
                    end
                    
                   if row(i) < controlrow45 && row(i) >= controlrow34
  
                        AAgent34 = find(seat==(row(i)*6-5));
                        BAgent34 = find(seat==(row(i)*6-4));
                        CAgent34 = find(seat==(row(i)*6-3));
                        DAgent34 = find(seat==(row(i)*6-2));
                        EAgent34 = find(seat==(row(i)*6-1));
                        FAgent34 = find(seat==(row(i)*6-0));
         
                        release_time(AAgent34) = start_group2 + (processT * randperm(num_group2,1)+ceil(limitcall*rand));
                        release_time(FAgent34) = start_group2 + (processT * randperm(num_group2,1)+ceil(limitcall*rand));
                        release_time(BAgent34) = start_group2 + (processT * randperm(num_group2,1)+ceil(limitcall*rand));
                        release_time(EAgent34) = start_group2 + (processT * randperm(num_group2,1)+ceil(limitcall*rand));
                        release_time(CAgent34) = start_group2 + (processT * randperm(num_group2,1)+ceil(limitcall*rand));
                        release_time(DAgent34) = start_group2 + (processT * randperm(num_group2,1)+ceil(limitcall*rand));
                        
                   end

                   if row(i) < controlrow34 && row(i) >= controlrow23
                       
                        AAgent23 = find(seat==(row(i)*6-5));
                        BAgent23 = find(seat==(row(i)*6-4));
                        CAgent23 = find(seat==(row(i)*6-3));
                        DAgent23 = find(seat==(row(i)*6-2));
                        EAgent23 = find(seat==(row(i)*6-1));
                        FAgent23 = find(seat==(row(i)*6-0));
                
                        release_time(AAgent23) = start_group2 + (processT * randperm(num_group2,1)+ceil(limitcall*rand));
                        release_time(FAgent23) = start_group2 + (processT * randperm(num_group2,1)+ceil(limitcall*rand));
                        release_time(BAgent23) = start_group2 + (processT * randperm(num_group2,1)+ceil(limitcall*rand));
                        release_time(EAgent23) = start_group2 + (processT * randperm(num_group2,1)+ceil(limitcall*rand));
                        release_time(CAgent23) = start_group2 + (processT * randperm(num_group2,1)+ceil(limitcall*rand));
                        release_time(DAgent23) = start_group2 + (processT * randperm(num_group2,1)+ceil(limitcall*rand));

                    end
                      
                   if row(i) < controlrow23 && row(i) >= controlrow12
                       
                       AAgent12 = find(seat==(row(i)*6-5));
                        BAgent12 = find(seat==(row(i)*6-4));
                        CAgent12 = find(seat==(row(i)*6-3));
                        DAgent12 = find(seat==(row(i)*6-2));
                        EAgent12 = find(seat==(row(i)*6-1));
                        FAgent12 = find(seat==(row(i)*6-0));
                  
                        release_time(AAgent12) = start_group2 + (processT * randperm(num_group2,1)+ceil(limitcall*rand));
                        release_time(FAgent12) = start_group2 + (processT * randperm(num_group2,1)+ceil(limitcall*rand));
                        release_time(BAgent12) = start_group2 + (processT * randperm(num_group2,1)+ceil(limitcall*rand));
                        release_time(EAgent12) = start_group2 + (processT * randperm(num_group2,1)+ceil(limitcall*rand));
                        release_time(CAgent12) = start_group3 + (processT * randperm(num_group3,1)+ceil(limitcall*rand));
                        release_time(DAgent12) = start_group3 + (processT * randperm(num_group3,1)+ceil(limitcall*rand));
 
                   end
                      
                   if row(i) < controlrow12
                        AAgent01 = find(seat==(row(i)*6-5));
                        BAgent01 = find(seat==(row(i)*6-4));
                        CAgent01 = find(seat==(row(i)*6-3));
                        DAgent01 = find(seat==(row(i)*6-2));
                        EAgent01 = find(seat==(row(i)*6-1));
                        FAgent01 = find(seat==(row(i)*6-0));

                        release_time(AAgent01) = start_group2 + (processT * randperm(num_group2,1)+ceil(limitcall*rand));
                        release_time(FAgent01) = start_group2 + (processT * randperm(num_group2,1)+ceil(limitcall*rand));
                        release_time(BAgent01) = start_group3 + (processT * randperm(num_group3,1)+ceil(limitcall*rand));
                        release_time(EAgent01) = start_group3 + (processT * randperm(num_group3,1)+ceil(limitcall*rand));
                        
                        if row(i) == ceil(controlrow12/4)
                            release_time(CAgent01) = start_group3 + (processT * (num_group3-1))+ceil(limitcall*rand);
                        else
                            release_time(CAgent01) = start_group3 + (processT * randperm(num_group3,1))+ceil(limitcall*rand);
                        end
                        if row(i) == ceil(3*controlrow12/4)
                            release_time(DAgent01) = start_group3 + (processT * (num_group3-2))+ceil(limitcall*rand);
                        else
                            release_time(DAgent01) = start_group3 + (processT * randperm(num_group3,1))+ceil(limitcall*rand);
                        end
                    end
        end    
end

%--- Main Loop Section (t>1) ---%
for t= 2:length(timelimit)
    for i = 1:crowdsize
       if ismember(i,M(:,1)) == 0
            Velocities(:,:,i,t) = Velocities(:,:,i,t-1) + 0.5.*Accelerations(:,:,i,t-1);
            Desiredspeed(i,t) = Desiredspeed(i,t-1);
            Maxspeed(i,t) = Desiredspeed(i,t)*1.3;
            if norm(Velocities(:,:,i,t)) > Maxspeed(i,t)
                arg = atan(Velocities(2,1,i,t)/Velocities(1,1,i,t));
                if Velocities(1,1,i,t) < 0
                    Velocities(:,:,i,t) = [-Maxspeed(i,t)*cos(arg);Maxspeed(i,t)*sin(arg)];
                else
                    Velocities(:,:,i,t) = [Maxspeed(i,t)*cos(arg);Maxspeed(i,t)*sin(arg)];
                end
            end
        end
    end
    for i = 1:crowdsize
            Coordinates(:,:,i,t) = Coordinates(:,:,i,t-1) + Velocities(:,:,i,t);
    end
    for i = 1:crowdsize
        if ismember(i,M) == 0
            for a = 1:crowdsize
            if ismember(a,M) == 0
                    if a ~= i
                        Ria = Coordinates(:,:,i,t) - Coordinates(:,:,a,t);
                        d = norm(Ria);
                        
                        if d <= PS
                            nRia = Ria./d;
                            nRiat = [-nRia(2);nRia(1)];
                            r = (diameter(i) + diameter(a))/2;
  
                               if dot(nRia, Velocities(:,:,i,t)./norm(Velocities(:,:,i,t))) < 0

                                            if norm(Velocities(:,:,i,t)) >= 0.6*Desiredspeed(i,t)
                                                Velocities(:,:,i,t) = Velocities(:,:,i,t).*d/PS;
                                            end

                                            if d <= r
                                                Forces(:,:,i,t) = Forces(:,:,i,t) + (A*exp((r-d)/B).*nRia + k*(r-d).*nRia + K*(r-d)*dot(Velocities(:,:,i,t),nRiat).*nRiat);
                                            else
                                                Forces(:,:,i,t) = Forces(:,:,i,t) + (A*exp((r-d)/B).*nRia);
                                            end

                               else            
                                            if d <= r
                                                Forces(:,:,i,t) = Forces(:,:,i,t) + 0.5.*(A*exp((r-d)/B).*nRia + k*(r-d).*nRia + K*(r-d)*dot(Velocities(:,:,i,t),nRiat).*nRiat);
                                            else
                                                Forces(:,:,i,t) = Forces(:,:,i,t) + 0.5.*(A*exp((r-d)/B).*nRia);
                                            end
                               end
                           
                        end
                    end
            end
            end
        end
    end
    for i = 1:crowdsize
    %--- Boundary for passenger ---%
        %--- Making sure passenger stays in pathway---%
        
            % Pathway (before entering the place)
            
                % Special border cases
                if Coordinates(2,1,i,t) >= pitch_total+front_space
                    Coordinates(2,1,i,t) = pitch_total+front_space-smallx;
                elseif Coordinates(2,1,i,t) <= pitch_total+1 && Coordinates(1,1,i,t) < 0
                    Coordinates(2,1,i,t) = pitch_total+1+smallx;
                end
                
                % Main Pathway cases
                
                if Coordinates(2,1,i,t) < pitch_total+front_space && Coordinates(2,1,i,t) > pitch_total+1
                    
                % Left boundary
                        if Coordinates(1,1,i,t) < spawnx
                            Coordinates(1,1,i,t) = spawnx+smallx;
                            
                % Right Boundary
                        elseif Coordinates(1,1,i,t) > widthX+width_walk
                            Coordinates(1,1,i,t) = widthX+width_walk-smallx;

                % People walk straight in the walkway
                        elseif Coordinates(1,1,i,t) < widthX && Coordinates(1,1,i,t) >= spawnx
                            %turn = [widthX+width_walk;pitch_total+1] - Coordinates(:,:,i,t);
                            %Rie = [widthX+width_walk;pitch_total+1] - Coordinates(:,:,i,t);
                            Desireddirections(:,:,i,t) = [1;0];
                            Forces(:,:,i,t) = Forces(:,:,i,t) + (1/T_i).*(Desiredspeed(i,t).*Desireddirections(:,:,i,t) - Velocities(:,:,i,t));
                            
                            % Walkway  wall force
                        	UiW =  Coordinates(:,:,i,t) - [Coordinates(1,1,i,t);pitch_total+front_space];
                            LiW =  Coordinates(:,:,i,t) - [Coordinates(1,1,i,t);pitch_total+1];
                            dU = norm(UiW);
                            dL = norm(LiW);
                            
                            r = diameter(i);  
                            nUiW = UiW./dU;
                            nUiWt = [-nUiW(2);nUiW(1)];
                            nLiW = LiW./dL;
                            nLiWt = [-nLiW(2);nLiW(1)];
                             
                             if dU <= r
                                Forces(:,:,i,t) = Forces(:,:,i,t) + UiW.*(U*exp(-dU/R)) + k*(r-dU).*nUiW - (K*(r-dU)*dot(Velocities(:,:,i,t),nUiWt)).*nUiWt;
                             else
                                Forces(:,:,i,t) = Forces(:,:,i,t) + UiW.*(U*exp(-dU/R));
                             end
                            
                              if dL <= r
                                Forces(:,:,i,t) = Forces(:,:,i,t) + LiW.*(U*exp(-dL/R)) + k*(r-dL).*nLiW - (K*(r-dL)*dot(Velocities(:,:,i,t),nLiWt)).*nLiWt; 
                              else
                                Forces(:,:,i,t) = Forces(:,:,i,t) + LiW.*(U*exp(-dL/R));
                             end
                                   
                % When there is a 90 degree turn, people need to start turning
                        elseif Coordinates(1,1,i,t) >= widthX && Coordinates(2,1,i,t) > pitch_total +1
                            
                            Desireddirections(:,:,i,t) = [0;-1];
                            Forces(:,:,i,t) = Forces(:,:,i,t) + (1/T_i).*(Desiredspeed(i,t).*Desireddirections(:,:,i,t) - Velocities(:,:,i,t));
                              
                         RiW =  Coordinates(:,:,i,t) - [widthX+width_walk;Coordinates(2,1,i,t)];
                         LiW =  Coordinates(:,:,i,t) - [widthX;Coordinates(2,1,i,t)];
                         UiW =  Coordinates(:,:,i,t) - [Coordinates(1,1,i,t);pitch_total+front_space];
                            dR = norm(RiW);
                            dL = norm(LiW);
                            dU = norm(UiW);
                            
                            r = diameter(i);  
                            nUiW = UiW./dU;
                            nUiWt = [-nUiW(2);nUiW(1)];
                            nLiW = LiW./dL;
                            nLiWt = [-nLiW(2);nLiW(1)];
                            nRiW = RiW./dR;
                            nRiWt = [-nRiW(2);nRiW(1)];
                            
                            if dU <= r
                                Forces(:,:,i,t) = Forces(:,:,i,t) + UiW.*(U*exp(-dU/R)) + k*(r-dU).*nUiW - (K*(r-dU)*dot(Velocities(:,:,i,t),nUiWt)).*nUiWt; 
                            else
                                Forces(:,:,i,t) = Forces(:,:,i,t) + UiW.*(U*exp(-dU/R));
                            end
                            
                              if dL <= r
                                Forces(:,:,i,t) = Forces(:,:,i,t) + LiW.*(U*exp(-dL/R)) + k*(r-dL).*nLiW - (K*(r-dL)*dot(Velocities(:,:,i,t),nLiWt)).*nLiWt;
                              else
                                Forces(:,:,i,t) = Forces(:,:,i,t) + LiW.*(U*exp(-dL/R));
                              end
                             
                              if dR <= r
                                Forces(:,:,i,t) = Forces(:,:,i,t) + RiW.*(U*exp(-dR/R)) + k*(r-dR).*nRiW - (K*(r-dR)*dot(Velocities(:,:,i,t),nRiWt)).*nRiWt;
                              else
                                Forces(:,:,i,t) = Forces(:,:,i,t) + RiW.*(U*exp(-dR/R));
                              end
                              
                                Rie = [widthX+width_walk/2;0] - Coordinates(:,:,i,t);
                                Desireddirections(:,:,i,t) = Rie./norm(Rie);
                                Forces(:,:,i,t) = Forces(:,:,i,t) + (1/T_i).*(Desiredspeed(i,t).*Desireddirections(:,:,i,t) - Velocities(:,:,i,t));
                        end
                end
              
                if Coordinates(2,1,i,t) <= pitch_total+1 && Coordinates(2,1,i,t) > pitch_total
                    
                    Desireddirections(:,:,i,t) = [0;-1];
                    Forces(:,:,i,t) = Forces(:,:,i,t) + (1/T_i).*(Desiredspeed(i,t).*Desireddirections(:,:,i,t) - Velocities(:,:,i,t));
                    
                     if Coordinates(1,1,i,t) < widthX
                           Coordinates(1,1,i,t) = widthX + smallx;
                            Velocities(1,1,i,t) = 0;
                     elseif Coordinates(1,1,i,t) > widthX + width_walk     
                           Coordinates(1,1,i,t) = widthX + width_walk - smallx;
                            Velocities(1,1,i,t) = 0;
                     end
          
                            RiW =  Coordinates(:,:,i,t) - [widthX+width_walk;Coordinates(2,1,i,t)];
                            LiW =  Coordinates(:,:,i,t) - [widthX;Coordinates(2,1,i,t)];
                            UiW =  Coordinates(:,:,i,t) - [Coordinates(1,1,i,t);pitch_total+front_space];
                            dR = norm(RiW);
                            dL = norm(LiW);
                            dU = norm(UiW);
                             
                            r = diameter(i);  
                            nUiW = UiW./dU;
                            nUiWt = [-nUiW(2);nUiW(1)];
                            nLiW = LiW./dL;
                            nLiWt = [-nLiW(2);nLiW(1)];
                            nRiW = RiW./dR;
                            nRiWt = [-nRiW(2);nRiW(1)];
                            
                            
                            if dU <= r
                                Forces(:,:,i,t) = Forces(:,:,i,t) + UiW.*(U*exp(-dU/R)) + k*(r-dU).*nUiW - (K*(r-dU)*dot(Velocities(:,:,i,t),nUiWt)).*nUiWt;
                            else
                                Forces(:,:,i,t) = Forces(:,:,i,t) + UiW.*(U*exp(-dU/R));
                             end
                            
                              if dL <= r
                                Forces(:,:,i,t) = Forces(:,:,i,t) + LiW.*(U*exp(-dL/R)) + k*(r-dL).*nLiW - (K*(r-dL)*dot(Velocities(:,:,i,t),nLiWt)).*nLiWt;
                              else
                                Forces(:,:,i,t) = Forces(:,:,i,t) + LiW.*(U*exp(-dL/R));
                              end
                             
                              if dR <= r
                               Forces(:,:,i,t) = Forces(:,:,i,t) + RiW.*(U*exp(-dR/R)) + k*(r-dR).*nRiW - (K*(r-dR)*dot(Velocities(:,:,i,t),nRiWt)).*nRiWt;
                              else
                               Forces(:,:,i,t) = Forces(:,:,i,t) + RiW.*(U*exp(-dR/R));
                              end
                                     
                end
     
            % On-board
       %--- Row pathways ---%

       
       for n=1:24
        if Coordinates(2,1,i,t) < pitch_total-pitchY*(n-1) && Coordinates(2,1,i,t) > pitch_total-n*pitchY && Coordinates(2,1,i,t) < pitch_total
            if Coordinates(1,1,i,t) < widthX || Coordinates(1,1,i,t) > width_total-width_walk
                    Coordinates(2,1,i,t) = pitch_total-(2*n-1)/2*pitchY;
                    %Velocities(2,1,i,t)=0;
            end
        end
       end

       
       
       if Coordinates(2,1,i,t) <= 0
           Coordinates(2,1,i,t) = smallx;
       end
       if Coordinates(1,1,i,t) < 0 && Coordinates(2,1,i,t) < pitch_total
           Coordinates(1,1,i,t) = smallx;
       end
       
                % Special border cases
                    if Coordinates(2,1,i,t) <= pitch_total && Coordinates(2,1,i,t) >= 0
                        
                            
                        if Coordinates(2,1,i,t) == pitch_total && Coordinates(1,1,i,t) > widthX+width_walk 
                             Coordinates(1,1,i,t) = widthX+width_walk-smallx;
                        end 
                        if Coordinates(2,1,i,t) == pitch_total && Coordinates(1,1,i,t) < widthX
                            Coordinates(1,1,i,t) = widthX+smallx;
                        end
                        
                        % Seating finding algo
                        if Coordinates(1,1,i,t) < 0
                                Coordinates(1,1,i,t) = smallx;
                                 Velocities(1,1,i,t) = 0;
                                 
                        elseif Coordinates(2,1,i,t) <= rowy(i) - pitchY/2
                            
                            if Coordinates(1,1,i,t) < widthX
                                Coordinates(1,1,i,t) = widthX+smallx;
                            end
                            if Coordinates(1,1,i,t) > widthX + width_walk
                                Coordinates(1,1,i,t) = widthX + width_walk-smallx;
                            end
                                    Desireddirections(:,:,i,t) = [0;1];
                                    Forces(:,:,i,t) = Forces(:,:,i,t) + (1/T_i).*(Desiredspeed(i,t).*Desireddirections(:,:,i,t) - Velocities(:,:,i,t));     
                        
                        elseif Coordinates(2,1,i,t) == rowy(i) + pitchY/2
                                if Coordinates(1,1,i,t) < widthX
                                    Coordinates(1,1,i,t) = widthX+smallx;
                                end
                                if Coordinates(1,1,i,t) > widthX + width_walk
                                    Coordinates(1,1,i,t) = widthX + width_walk-smallx;
                                end
                                    Desireddirections(:,:,i,t) = [0;-1];
                                    Forces(:,:,i,t) = Forces(:,:,i,t) + (1/T_i).*(Desiredspeed(i,t).*Desireddirections(:,:,i,t) - Velocities(:,:,i,t));
                                    
                        elseif   Coordinates(2,1,i,t) < rowy(i) + pitchY/2 && Coordinates(2,1,i,t) > rowy(i) - pitchY/2 && Coordinates(1,1,i,t) <= 0   
                                Coordinates(:,1,i,t) = [widthX/6;rowy(i)];
                                 
                        elseif   Coordinates(2,1,i,t) < rowy(i) + pitchY/2 && Coordinates(2,1,i,t) > rowy(i) - pitchY/2 && Coordinates(1,1,i,t) > 0
                            
                                Desiredspeed(i,t) = Desiredspeed(i,1)/3.3;

                                %Column A algo
                                if mod(seat(i),6) == 1
                                            if Coordinates(1,1,i,t) < 0
                                                 Coordinates(1,1,i,t) = smallx;
                                             elseif Coordinates(1,1,i,t) >  2*widthX + width_walk;
                                                 Coordinates(1,1,i,t) = widthX + width_walk - smallx;
                                             end
                                            %Rie = [widthX/3;rowy(i)] - Coordinates(:,:,i,t);
                                            Desireddirections(:,:,i,t) = [-1;0];
                                            Forces(:,:,i,t) = Forces(:,:,i,t) + (1/T_i).*(Desiredspeed(i,t).*Desireddirections(:,:,i,t) - Velocities(:,:,i,t));

                                            AAgent = find(seat==(row(i)*6-5));
                                            BAgent = find(seat==(row(i)*6-4));
                                            CAgent = find(seat==(row(i)*6-3));
                                            
                                            if luggage(AAgent) ==1
                                                if t <= stop_lug(AAgent) && stop_lug(AAgent) ~= 0
                                                    Velocities(:,:,AAgent,t) = 0;
                                                    Coordinates(:,1,i,t) = [(widthX+smallx);rowy(i)+ 0.1*(pitchY/2)*5/6];
                                                    waitlug(AAgent) = 1;
                                                end
                                                
                                                if start_lug(AAgent) == 0 && stop_lug(AAgent) == 0
                                                start_lug(AAgent) = t;
                                                stop_lug(AAgent) = t+luggage_time;
                                                end
                                            end
                                            
                                         if Coordinates(1,1,AAgent,t) <= (10*widthX/36) && Coordinates(2,1,AAgent,t) < pitch_total 
                                              if inseat(AAgent) == 0 && t > release_time(AAgent)
                                                  inseat(AAgent) = t;
                                              end
                                                Velocities(:,:,i,t) = 0;
                                                Desireddirections(:,:,i,t) = [0;0];
                                                Coordinates(:,1,i,t) = [(10*widthX/36);rowy(i)];
                                                blocknmoving(AAgent) = 0;
                                          end
  
                                        if Coordinates(1,1,i,t) > (10*widthX/36)
                                            
                                            if blocknmoving(AAgent) == 1 && blocknmoving(BAgent) == 1
                                                Coordinates(:,1,BAgent,t) = [widthX+width_walk-smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                            end
                                            if blocknmoving(AAgent) == 1 && blocknmoving(CAgent) == 1
                                                Coordinates(:,1,CAgent,t) = [widthX+width_walk-smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                            end
                                            if blocknmoving(BAgent) == 1 && blocknmoving(CAgent) == 1
                                                Coordinates(:,1,CAgent,t) = [widthX+width_walk-smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                            end
                                            if blocknmoving(AAgent) == 1 && blocknmoving(BAgent) == 1  && blocknmoving(CAgent) == 1
                                                Coordinates(:,1,BAgent,t) = [widthX+width_walk-2*smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                                Coordinates(:,1,CAgent,t) = [widthX+width_walk-smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                            end
                                              
                                            % If Column C is occupied
                                            if inseat(CAgent) ~= 0 && inseat(AAgent) == 0 && inseat(BAgent) == 0
                                                if or(Coordinates(2,1,AAgent,t) > rowy(i) - pitchY/2 && Coordinates(2,1,AAgent,t) < rowy(i) + pitchY/2, Coordinates(2,1,BAgent,t) > rowy(i) - pitchY/2 && Coordinates(2,1,BAgent,t) < rowy(i) + pitchY/2)
                                                   Coordinates(:,1,CAgent,t) = [widthX+width_walk-smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                                   blocknmoving(CAgent) = 1;
                                                   inseat(CAgent) = 0;
                                                   seat1_interference =  seat1_interference+1;
                                                end
                                            %If Column B is occupied
                                            elseif inseat(BAgent) ~= 0 && inseat(AAgent) == 0 && inseat(CAgent) == 0
                                                if Coordinates(2,1,AAgent,t) > rowy(i) - pitchY/2 && Coordinates(2,1,AAgent,t) < rowy(i) + pitchY/2
                                                   Coordinates(:,1,BAgent,t) = [widthX+width_walk-smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                                   blocknmoving(BAgent) = 1;
                                                   inseat(BAgent) = 0;
                                                   seat1_interference =  seat1_interference+1;
                                                end
                                            %If Column B and  C are occupied
                                            elseif inseat(BAgent) ~= 0 && inseat(CAgent) ~= 0 && inseat(AAgent) == 0 && Coordinates(2,1,AAgent,t) > rowy(i) - pitchY/2 && Coordinates(2,1,AAgent,t) < rowy(i) + pitchY/2
                                               Coordinates(:,1,BAgent,t) = [widthX+width_walk-2*smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                               Coordinates(:,1,CAgent,t) = [widthX+width_walk-smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                               blocknmoving(BAgent) = 1;
                                               blocknmoving(CAgent) = 1;
                                               inseat(BAgent) = 0;
                                               inseat(CAgent) = 0;
                                               seat2_interference =  seat2_interference+1;
                                            %If Column A and C are occupied
                                            elseif inseat(AAgent) ~= 0 && inseat(CAgent) ~= 0 && inseat(BAgent) == 0 && Coordinates(2,1,BAgent,t) > rowy(i) - pitchY/2 && Coordinates(2,1,BAgent,t) < rowy(i) + pitchY/2
                                               Coordinates(:,1,CAgent,t) = [widthX+width_walk-smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                               blocknmoving(CAgent) = 1;
                                               blocknmoving(AAgent) = 0;
                                               inseat(CAgent) = 0;
                                               seat1_interference =  seat1_interference+1;
                                            end
                                        end
                                        
                                    
                                 %Column B algo                  
                                elseif mod(seat(i),6) == 2
                                    if Coordinates(1,1,i,t) < 0
                                        Coordinates(1,1,i,t) = smallx;
                                     elseif Coordinates(1,1,i,t) >  2*widthX + width_walk;
                                         Coordinates(1,1,i,t) = widthX + width_walk - smallx;
                                     end
                                    %Rie = [widthX/2;rowy(i)] - Coordinates(:,:,i,t);
                                    Desireddirections(:,:,i,t) = [-1;0];
                                    Forces(:,:,i,t) = Forces(:,:,i,t) + (1/T_i).*(Desiredspeed(i,t).*Desireddirections(:,:,i,t) - Velocities(:,:,i,t));
                                    
                                    
                                            AAgent = find(seat==(row(i)*6-5));
                                            BAgent = find(seat==(row(i)*6-4));
                                            CAgent = find(seat==(row(i)*6-3));
                                            
                                            if luggage(BAgent) ==1
                                                if t <= stop_lug(BAgent) && stop_lug(BAgent) ~= 0
                                                    Velocities(:,:,BAgent,t) = 0;
                                                    Coordinates(:,1,i,t) = [(widthX+smallx);rowy(i)+ 0.1*(pitchY/2)*3/6];
                                                end
                                                
                                                if start_lug(BAgent) == 0 && stop_lug(BAgent) == 0
                                                start_lug(BAgent) = t;
                                                stop_lug(BAgent) = t+luggage_time;
                                                end
                                            end
                                            
                                        if Coordinates(1,1,BAgent,t) <= (20*widthX/36) && Coordinates(2,1,BAgent,t) < pitch_total
                                            if inseat(BAgent) == 0 && t > release_time(BAgent)
                                                  inseat(BAgent) = t;
                                            end 
                                            Velocities(:,:,i,t)=0;
                                            Coordinates(:,1,i,t) = [(20*widthX/36);rowy(i)];
                                            blocknmoving(BAgent) = 0;
                                            Forces(:,:,i,t) = 0;
                                        end
          
                                        if Coordinates(1,1,i,t) > (20*widthX/36)
                                            
                                            if blocknmoving(AAgent) == 1 && blocknmoving(BAgent) == 1
                                                Coordinates(:,1,BAgent,t) = [widthX+width_walk-smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                            end
                                            if blocknmoving(AAgent) == 1 && blocknmoving(CAgent) == 1
                                                Coordinates(:,1,CAgent,t) = [widthX+width_walk-smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                            end
                                            if blocknmoving(BAgent) == 1 && blocknmoving(CAgent) == 1
                                                Coordinates(:,1,CAgent,t) = [widthX+width_walk-smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                            end
                                            if blocknmoving(AAgent) == 1 && blocknmoving(BAgent) == 1  && blocknmoving(CAgent) == 1
                                                Coordinates(:,1,BAgent,t) = [widthX+width_walk-2*smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                                Coordinates(:,1,CAgent,t) = [widthX+width_walk-smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                            end
                                              
                                            % If Column C is occupied
                                            if inseat(CAgent) ~= 0 && inseat(AAgent) == 0 && inseat(BAgent) == 0
                                                if or(Coordinates(2,1,AAgent,t) > rowy(i) - pitchY/2 && Coordinates(2,1,AAgent,t) < rowy(i) + pitchY/2, Coordinates(2,1,BAgent,t) > rowy(i) - pitchY/2 && Coordinates(2,1,BAgent,t) < rowy(i) + pitchY/2)
                                                   Coordinates(:,1,CAgent,t) = [widthX+width_walk-smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                                   blocknmoving(CAgent) = 1;
                                                   inseat(CAgent) = 0;
                                                   seat1_interference =  seat1_interference+1;
                                                end
                                            %If Column B is occupied
                                            elseif inseat(BAgent) ~= 0 && inseat(AAgent) == 0 && inseat(CAgent) == 0
                                                if Coordinates(2,1,AAgent,t) > rowy(i) - pitchY/2 && Coordinates(2,1,AAgent,t) < rowy(i) + pitchY/2
                                                   Coordinates(:,1,BAgent,t) = [widthX+width_walk-smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                                   blocknmoving(BAgent) = 1;
                                                   inseat(BAgent) = 0;
                                                   seat1_interference =  seat1_interference+1;
                                                end
                                            %If Column B and  C are occupied
                                            elseif inseat(BAgent) ~= 0 && inseat(CAgent) ~= 0 && inseat(AAgent) == 0 && Coordinates(2,1,AAgent,t) > rowy(i) - pitchY/2 && Coordinates(2,1,AAgent,t) < rowy(i) + pitchY/2
                                               Coordinates(:,1,BAgent,t) = [widthX+width_walk-2*smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                               Coordinates(:,1,CAgent,t) = [widthX+width_walk-smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                               blocknmoving(BAgent) = 1;
                                               blocknmoving(CAgent) = 1;
                                               inseat(BAgent) = 0;
                                               inseat(CAgent) = 0;
                                               seat2_interference =  seat2_interference+1;
                                            %If Column A and C are occupied
                                            elseif inseat(AAgent) ~= 0 && inseat(CAgent) ~= 0 && inseat(BAgent) == 0 && Coordinates(2,1,BAgent,t) > rowy(i) - pitchY/2 && Coordinates(2,1,BAgent,t) < rowy(i) + pitchY/2
                                               Coordinates(:,1,CAgent,t) = [widthX+width_walk-smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                               blocknmoving(CAgent) = 1;
                                               blocknmoving(AAgent) = 0;
                                               inseat(CAgent) = 0;
                                               seat1_interference =  seat1_interference+1;
                                            end
                                        end

                                    
                                elseif mod(seat(i),6) == 3
                                    
                                    if Coordinates(1,1,i,t) < 0
                                        Coordinates(1,1,i,t) = smallx;
                                     elseif Coordinates(1,1,i,t) >  2*widthX + width_walk;
                                         Coordinates(1,1,i,t) = widthX + width_walk - smallx;
                                    end
                                    %Rie = [5/6*widthX;rowy(i)] - Coordinates(:,:,i,t);
                                    Desireddirections(:,:,i,t) = [-1;0];
                                    Forces(:,:,i,t) = Forces(:,:,i,t) + (1/T_i).*(Desiredspeed(i,t).*Desireddirections(:,:,i,t) - Velocities(:,:,i,t));        
                                    
                                            AAgent = find(seat==(row(i)*6-5));
                                            BAgent = find(seat==(row(i)*6-4));
                                            CAgent = find(seat==(row(i)*6-3));
                                            
                                            if luggage(CAgent) ==1
                                                if t <= stop_lug(CAgent) && stop_lug(CAgent) ~= 0
                                                    Velocities(:,:,CAgent,t) = 0;
                                                    Coordinates(:,1,i,t) = [(widthX+smallx);rowy(i)+ 0.1*(pitchY/2)*1/6];
                                                end
                                                
                                                if start_lug(CAgent) == 0 && stop_lug(CAgent) == 0
                                                start_lug(CAgent) = t;
                                                stop_lug(CAgent) = t+luggage_time;
                                                end
                                            end

                                         if Coordinates(1,1,CAgent,t) <= (30*widthX/36) && Coordinates(2,1,CAgent,t) < pitch_total
                                            if inseat(CAgent) == 0 && t > release_time(CAgent)
                                                  inseat(CAgent) = t;
                                            end 
                                           Velocities(:,:,i,t)=0;
                                           Coordinates(:,1,i,t) = [(30*widthX/36);rowy(i)];
                                           blocknmoving(CAgent) = 0;
                                           Forces(:,:,i,t) = 0;     
                                         end
                                                         
                                        if Coordinates(1,1,i,t) > (30*widthX/36)
                                            
                                            if blocknmoving(AAgent) == 1 && blocknmoving(BAgent) == 1
                                                Coordinates(:,1,BAgent,t) = [widthX+width_walk-smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                            end
                                            if blocknmoving(AAgent) == 1 && blocknmoving(CAgent) == 1
                                                Coordinates(:,1,CAgent,t) = [widthX+width_walk-smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                            end
                                            if blocknmoving(BAgent) == 1 && blocknmoving(CAgent) == 1
                                                Coordinates(:,1,CAgent,t) = [widthX+width_walk-smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                            end
                                            if blocknmoving(AAgent) == 1 && blocknmoving(BAgent) == 1  && blocknmoving(CAgent) == 1
                                                Coordinates(:,1,BAgent,t) = [widthX+width_walk-2*smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                                Coordinates(:,1,CAgent,t) = [widthX+width_walk-smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                            end
                                              
                                            % If Column C is occupied
                                            if inseat(CAgent) ~= 0 && inseat(AAgent) == 0 && inseat(BAgent) == 0
                                                if or(Coordinates(2,1,AAgent,t) > rowy(i) - pitchY/2 && Coordinates(2,1,AAgent,t) < rowy(i) + pitchY/2, Coordinates(2,1,BAgent,t) > rowy(i) - pitchY/2 && Coordinates(2,1,BAgent,t) < rowy(i) + pitchY/2)
                                                   Coordinates(:,1,CAgent,t) = [widthX+width_walk-smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                                   blocknmoving(CAgent) = 1;
                                                   inseat(CAgent) = 0;
                                                   seat1_interference =  seat1_interference+1;
                                                end
                                            %If Column B is occupied
                                            elseif inseat(BAgent) ~= 0 && inseat(AAgent) == 0 && inseat(CAgent) == 0
                                                if Coordinates(2,1,AAgent,t) > rowy(i) - pitchY/2 && Coordinates(2,1,AAgent,t) < rowy(i) + pitchY/2
                                                   Coordinates(:,1,BAgent,t) = [widthX+width_walk-smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                                   blocknmoving(BAgent) = 1;
                                                   inseat(BAgent) = 0;
                                                   seat1_interference =  seat1_interference+1;
                                                end
                                            %If Column B and  C are occupied
                                            elseif inseat(BAgent) ~= 0 && inseat(CAgent) ~= 0 && inseat(AAgent) == 0 && Coordinates(2,1,AAgent,t) > rowy(i) - pitchY/2 && Coordinates(2,1,AAgent,t) < rowy(i) + pitchY/2
                                               Coordinates(:,1,BAgent,t) = [widthX+width_walk-2*smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                               Coordinates(:,1,CAgent,t) = [widthX+width_walk-smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                               blocknmoving(BAgent) = 1;
                                               blocknmoving(CAgent) = 1;
                                               inseat(BAgent) = 0;
                                               inseat(CAgent) = 0;
                                               seat2_interference =  seat2_interference+1;
                                            %If Column A and C are occupied
                                            elseif inseat(AAgent) ~= 0 && inseat(CAgent) ~= 0 && inseat(BAgent) == 0 && Coordinates(2,1,BAgent,t) > rowy(i) - pitchY/2 && Coordinates(2,1,BAgent,t) < rowy(i) + pitchY/2
                                               Coordinates(:,1,CAgent,t) = [widthX+width_walk-smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                               blocknmoving(CAgent) = 1;
                                               blocknmoving(AAgent) = 0;
                                               inseat(CAgent) = 0;
                                               seat1_interference =  seat1_interference+1;
                                            end
                                        end
                                        
                                        
                                elseif mod(seat(i),6) == 4
                                    if Coordinates(1,1,i,t) < 0
                                        Coordinates(1,1,i,t) = smallx;
                                     elseif Coordinates(1,1,i,t) >  2*widthX + width_walk;
                                         Coordinates(1,1,i,t) = 2*widthX + width_walk - smallx;
                                    end
                                    Desireddirections(:,:,i,t) = [1;0];
                                    Forces(:,:,i,t) = Forces(:,:,i,t) + (1/T_i).*(Desiredspeed(i,t).*Desireddirections(:,:,i,t) - Velocities(:,:,i,t));  
                                    
                                            DAgent = find(seat==(row(i)*6-2));
                                            EAgent = find(seat==(row(i)*6-1));
                                            FAgent = find(seat==(row(i)*6-0));

                                             if luggage(DAgent) ==1
                                                if t <= stop_lug(DAgent) && stop_lug(DAgent) ~= 0
                                                    Velocities(:,:,DAgent,t) = 0;
                                                    Coordinates(:,1,i,t) = [(widthX+width_walk-smallx);rowy(i)- 0.1*(pitchY/2)*1/6];
                                                end
                                                if start_lug(DAgent) == 0 && stop_lug(DAgent) == 0
                                                start_lug(DAgent) = t;
                                                stop_lug(DAgent) = t+luggage_time;
                                                end
                                             end
                                             
                                         if Coordinates(1,1,DAgent,t) >=  (widthX+width_walk+(6*widthX/36)) && Coordinates(2,1,DAgent,t) < pitch_total
                                            if inseat(DAgent) == 0 && t > release_time(DAgent)
                                                  inseat(DAgent) = t;
                                            end 
                                           Velocities(:,:,DAgent,t)=0;
                                           Coordinates(:,:,DAgent,t) = [widthX+width_walk+(6*widthX/36);rowy(i)];
                                           blocknmoving(DAgent) = 0;
                                         end

                                        if Coordinates(1,1,DAgent,t) < (widthX+width_walk+(6*widthX/36))

                                            if blocknmoving(FAgent) == 1 && blocknmoving(EAgent) == 1
                                                Coordinates(:,1,EAgent,t) = [widthX+smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                            end
                                            if blocknmoving(FAgent) == 1 && blocknmoving(DAgent) == 1
                                                Coordinates(:,1,DAgent,t) = [widthX+smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                            end
                                            if blocknmoving(EAgent) == 1 && blocknmoving(DAgent) == 1
                                                Coordinates(:,1,DAgent,t) = [widthX+smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                            end
                                            if blocknmoving(FAgent) == 1 && blocknmoving(EAgent) == 1  && blocknmoving(DAgent) == 1
                                                Coordinates(:,1,EAgent,t) = [widthX+2*smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                                Coordinates(:,1,DAgent,t) = [widthX+smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                            end
                                            
                                            if inseat(FAgent) == 0 && Coordinates(2,1,FAgent,t) > rowy(i) - pitchY/2 && Coordinates(2,1,FAgent,t) < rowy(i) + pitchY/2
                                                if inseat(EAgent) ~= 0
                                                   Coordinates(:,1,EAgent,t) = [widthX+smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                                   blocknmoving(EAgent) = 1;
                                                   inseat(EAgent) = 0;
                                                end
                                                if inseat(DAgent) ~= 0
                                                   Coordinates(:,1,DAgent,t) = [widthX+smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                                   blocknmoving(DAgent) = 1;
                                                   inseat(DAgent) = 0;
                                                end
                                            end
                                            
                                            % If Column D is occupied
                                            if inseat(DAgent) ~= 0 && inseat(FAgent) == 0 && inseat(EAgent) == 0
                                                if or(Coordinates(2,1,FAgent,t) > rowy(i) - pitchY/2 && Coordinates(2,1,FAgent,t) < rowy(i) + pitchY/2, Coordinates(2,1,EAgent,t) > rowy(i) - pitchY/2 && Coordinates(2,1,EAgent,t) < rowy(i) + pitchY/2)
                                                   Coordinates(:,1,DAgent,t) = [widthX+smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                                   blocknmoving(DAgent) = 1;
                                                   inseat(DAgent) = 0;
                                                   seat1_interference =  seat1_interference+1;
                                                end
                                            %If Column E is occupied
                                            elseif inseat(EAgent) ~= 0 && inseat(FAgent) == 0 && inseat(DAgent) == 0
                                                if Coordinates(2,1,FAgent,t) > rowy(i) - pitchY/2 && Coordinates(2,1,FAgent,t) < rowy(i) + pitchY/2
                                                   Coordinates(:,1,EAgent,t) = [widthX+smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                                   blocknmoving(EAgent) = 1;
                                                   inseat(EAgent) = 0;
                                                   seat1_interference =  seat1_interference+1;
                                                end
                                            %If Column E and  D are occupied
                                            elseif inseat(EAgent) ~= 0 && inseat(DAgent) ~= 0 && inseat(FAgent) == 0 && Coordinates(2,1,FAgent,t) > rowy(i) - pitchY/2 && Coordinates(2,1,FAgent,t) < rowy(i) + pitchY/2
                                               Coordinates(:,1,EAgent,t) = [widthX+2*smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                               Coordinates(:,1,DAgent,t) = [widthX+smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                               blocknmoving(EAgent) = 1;
                                               blocknmoving(DAgent) = 1;
                                               inseat(EAgent) = 0;
                                               inseat(DAgent) = 0;
                                               seat2_interference =  seat2_interference+1;
                                            %If Column F and D are occupied
                                            elseif inseat(FAgent) ~= 0 && inseat(DAgent) ~= 0 && inseat(EAgent) == 0 && Coordinates(2,1,EAgent,t) > rowy(i) - pitchY/2 && Coordinates(2,1,EAgent,t) < rowy(i) + pitchY/2
                                               Coordinates(:,1,DAgent,t) = [widthX+smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                               blocknmoving(DAgent) = 1;
                                               blocknmoving(FAgent) = 0;
                                               inseat(DAgent) = 0;
                                               seat1_interference =  seat1_interference+1;
                                            end
                                       end
                                    
                                elseif mod(seat(i),6) == 5
                                    if Coordinates(1,1,i,t) < 0
                                        Coordinates(1,1,i,t) = smallx;
                                     elseif Coordinates(1,1,i,t) >  2*widthX + width_walk;
                                         Coordinates(1,1,i,t) = 2*widthX + width_walk - smallx;
                                     end
                                    %Rie = [widthX+width_walk+(16*widthX/36);rowy(i)] - Coordinates(:,:,i,t);
                                    Desireddirections(:,:,i,t) = [1;0];
                                    Forces(:,:,i,t) = Forces(:,:,i,t) + (1/T_i).*(Desiredspeed(i,t).*Desireddirections(:,:,i,t) - Velocities(:,:,i,t));   
                                    
                                            DAgent = find(seat==(row(i)*6-2));
                                            EAgent = find(seat==(row(i)*6-1));
                                            FAgent = find(seat==(row(i)*6-0));
                                            
                                             if luggage(EAgent) ==1
                                                if t <= stop_lug(EAgent) && stop_lug(EAgent) ~= 0
                                                    Velocities(:,:,EAgent,t) = 0;
                                                    Coordinates(:,1,i,t) = [(widthX+width_walk-smallx);rowy(i)- 0.1*(pitchY/2)*3/6];
                                                end
                                                if start_lug(EAgent) == 0 && stop_lug(EAgent) == 0
                                                start_lug(EAgent) = t;
                                                stop_lug(EAgent) = t+luggage_time;
                                                end
                                             end
       
                                          if Coordinates(1,1,EAgent,t) >=  (widthX+width_walk+(16*widthX/36)) && Coordinates(2,1,EAgent,t) < pitch_total
                                            if inseat(EAgent) == 0 && t > release_time(EAgent)
                                                  inseat(EAgent) = t;
                                            end 
                                           Velocities(:,:,EAgent,t)=0;
                                           Coordinates(:,:,EAgent,t) = [widthX+width_walk+(16*widthX/36);rowy(i)];
                                           blocknmoving(EAgent) = 0;
                                         end

                                       if Coordinates(1,1,EAgent,t) < (widthX+width_walk+(16*widthX/36))
                                           
                                          if blocknmoving(FAgent) == 1 && blocknmoving(EAgent) == 1
                                                Coordinates(:,1,EAgent,t) = [widthX+smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                          end
                                            if blocknmoving(FAgent) == 1 && blocknmoving(DAgent) == 1
                                                Coordinates(:,1,DAgent,t) = [widthX+smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                            end
                                            if blocknmoving(EAgent) == 1 && blocknmoving(DAgent) == 1
                                                Coordinates(:,1,DAgent,t) = [widthX+smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                            end
                                            if blocknmoving(FAgent) == 1 && blocknmoving(EAgent) == 1  && blocknmoving(DAgent) == 1
                                                Coordinates(:,1,EAgent,t) = [widthX+2*smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                                Coordinates(:,1,DAgent,t) = [widthX+smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                            end
                                            
                                            if inseat(FAgent) == 0 && Coordinates(2,1,FAgent,t) > rowy(i) - pitchY/2 && Coordinates(2,1,FAgent,t) < rowy(i) + pitchY/2
                                                if inseat(EAgent) ~= 0
                                                   Coordinates(:,1,EAgent,t) = [widthX+smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                                   blocknmoving(EAgent) = 1;
                                                   inseat(EAgent) = 0;
                                                end
                                                if inseat(DAgent) ~= 0
                                                   Coordinates(:,1,DAgent,t) = [widthX+smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                                   blocknmoving(DAgent) = 1;
                                                   inseat(DAgent) = 0;
                                                end
                                            end
                                            
                                            % If Column D is occupied
                                            if inseat(DAgent) ~= 0 && inseat(FAgent) == 0 && inseat(EAgent) == 0
                                                if or(Coordinates(2,1,FAgent,t) > rowy(i) - pitchY/2 && Coordinates(2,1,FAgent,t) < rowy(i) + pitchY/2, Coordinates(2,1,EAgent,t) > rowy(i) - pitchY/2 && Coordinates(2,1,EAgent,t) < rowy(i) + pitchY/2)
                                                   Coordinates(:,1,DAgent,t) = [widthX+smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                                   blocknmoving(DAgent) = 1;
                                                   inseat(DAgent) = 0;
                                                   seat1_interference =  seat1_interference+1;
                                                end
                                            %If Column E is occupied
                                            elseif inseat(EAgent) ~= 0 && inseat(FAgent) == 0 && inseat(DAgent) == 0
                                                if Coordinates(2,1,FAgent,t) > rowy(i) - pitchY/2 && Coordinates(2,1,FAgent,t) < rowy(i) + pitchY/2
                                                   Coordinates(:,1,EAgent,t) = [widthX+smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                                   blocknmoving(EAgent) = 1;
                                                   inseat(EAgent) = 0;
                                                   seat1_interference =  seat1_interference+1;
                                                end
                                            %If Column E and  D are occupied
                                            elseif inseat(EAgent) ~= 0 && inseat(DAgent) ~= 0 && inseat(FAgent) == 0 && Coordinates(2,1,FAgent,t) > rowy(i) - pitchY/2 && Coordinates(2,1,FAgent,t) < rowy(i) + pitchY/2
                                               Coordinates(:,1,EAgent,t) = [widthX+2*smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                               Coordinates(:,1,DAgent,t) = [widthX+smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                               blocknmoving(EAgent) = 1;
                                               blocknmoving(DAgent) = 1;
                                               inseat(EAgent) = 0;
                                               inseat(DAgent) = 0;
                                               seat2_interference =  seat2_interference+1;
                                            %If Column F and D are occupied
                                            elseif inseat(FAgent) ~= 0 && inseat(DAgent) ~= 0 && inseat(EAgent) == 0 && Coordinates(2,1,EAgent,t) > rowy(i) - pitchY/2 && Coordinates(2,1,EAgent,t) < rowy(i) + pitchY/2
                                               Coordinates(:,1,DAgent,t) = [widthX+smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                               blocknmoving(DAgent) = 1;
                                               blocknmoving(FAgent) = 0;
                                               inseat(DAgent) = 0;
                                               seat1_interference =  seat1_interference+1;
                                            end
                                       end
              
                                elseif mod(seat(i),6) == 0
                                    
                                    if Coordinates(1,1,i,t) < 0
                                        Coordinates(1,1,i,t) = smallx;
                                     elseif Coordinates(1,1,i,t) >  2*widthX + width_walk;
                                         Coordinates(1,1,i,t) = 2*widthX + width_walk - smallx;
                                     end
                                    %Rie = [widthX+width_walk+(26*widthX/36);rowy(i)] - Coordinates(:,:,i,t);
                                    Desireddirections(:,:,i,t) = [1;0];
                                    Forces(:,:,i,t) = Forces(:,:,i,t) + (1/T_i).*(Desiredspeed(i,t).*Desireddirections(:,:,i,t) - Velocities(:,:,i,t));   
                                    
                                            DAgent = find(seat==(row(i)*6-2));
                                            EAgent = find(seat==(row(i)*6-1));
                                            FAgent = find(seat==(row(i)*6-0));
                                            
                                             if luggage(FAgent) ==1
                                                if t <= stop_lug(FAgent) && stop_lug(FAgent) ~= 0
                                                    Velocities(:,:,FAgent,t) = 0;
                                                    Coordinates(:,1,i,t) = [(widthX+width_walk-smallx);rowy(i)- 0.1*(pitchY/2)*5/6];
                                                end
                                                if start_lug(FAgent) == 0 && stop_lug(FAgent) == 0
                                                start_lug(FAgent) = t;
                                                stop_lug(FAgent) = t+luggage_time;
                                                end
                                             end
        
                                           if Coordinates(1,1,FAgent,t) >=  (widthX+width_walk+(26*widthX/36)) && Coordinates(2,1,FAgent,t) < pitch_total
                                            if inseat(FAgent) == 0 && t > release_time(FAgent)
                                                  inseat(FAgent) = t;
                                            end 
                                           Velocities(:,:,FAgent,t)=0;
                                           Coordinates(:,:,FAgent,t) = [widthX+width_walk+(26*widthX/36);rowy(i)];
                                           blocknmoving(FAgent) = 0;
                                           end
              
                                       if Coordinates(1,1,FAgent,t) < (widthX+width_walk+(26*widthX/36))
                                            
                                           if blocknmoving(FAgent) == 1 && blocknmoving(EAgent) == 1
                                                Coordinates(:,1,EAgent,t) = [widthX+smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                            end
                                            if blocknmoving(FAgent) == 1 && blocknmoving(DAgent) == 1
                                                Coordinates(:,1,DAgent,t) = [widthX+smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                            end
                                            if blocknmoving(EAgent) == 1 && blocknmoving(DAgent) == 1
                                                Coordinates(:,1,DAgent,t) = [widthX+smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                            end
                                            if blocknmoving(FAgent) == 1 && blocknmoving(EAgent) == 1  && blocknmoving(DAgent) == 1
                                                Coordinates(:,1,EAgent,t) = [widthX+2*smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                                Coordinates(:,1,DAgent,t) = [widthX+smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                            end
                                            
                                            if inseat(FAgent) == 0 && Coordinates(2,1,FAgent,t) > rowy(i) - pitchY/2 && Coordinates(2,1,FAgent,t) < rowy(i) + pitchY/2
                                                if inseat(EAgent) ~= 0
                                                   Coordinates(:,1,EAgent,t) = [widthX+smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                                   blocknmoving(EAgent) = 1;
                                                   inseat(EAgent) = 0;
                                                end
                                                if inseat(DAgent) ~= 0
                                                   Coordinates(:,1,DAgent,t) = [widthX+smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                                   blocknmoving(DAgent) = 1;
                                                   inseat(DAgent) = 0;
                                                end
                                            end
                                            
                                            % If Column D is occupied
                                            if inseat(DAgent) ~= 0 && inseat(FAgent) == 0 && inseat(EAgent) == 0
                                                if or(Coordinates(2,1,FAgent,t) > rowy(i) - pitchY/2 && Coordinates(2,1,FAgent,t) < rowy(i) + pitchY/2, Coordinates(2,1,EAgent,t) > rowy(i) - pitchY/2 && Coordinates(2,1,EAgent,t) < rowy(i) + pitchY/2)
                                                   Coordinates(:,1,DAgent,t) = [widthX+smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                                   blocknmoving(DAgent) = 1;
                                                   inseat(DAgent) = 0;
                                                   seat1_interference =  seat1_interference+1;
                                                end
                                            %If Column E is occupied
                                            elseif inseat(EAgent) ~= 0 && inseat(FAgent) == 0 && inseat(DAgent) == 0
                                                if Coordinates(2,1,FAgent,t) > rowy(i) - pitchY/2 && Coordinates(2,1,FAgent,t) < rowy(i) + pitchY/2
                                                   Coordinates(:,1,EAgent,t) = [widthX+smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                                   blocknmoving(EAgent) = 1;
                                                   inseat(EAgent) = 0;
                                                   seat1_interference =  seat1_interference+1;
                                                end
                                            %If Column E and  D are occupied
                                            elseif inseat(EAgent) ~= 0 && inseat(DAgent) ~= 0 && inseat(FAgent) == 0 && Coordinates(2,1,FAgent,t) > rowy(i) - pitchY/2 && Coordinates(2,1,FAgent,t) < rowy(i) + pitchY/2
                                               Coordinates(:,1,EAgent,t) = [widthX+2*smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                               Coordinates(:,1,DAgent,t) = [widthX+smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                               blocknmoving(EAgent) = 1;
                                               blocknmoving(DAgent) = 1;
                                               inseat(EAgent) = 0;
                                               inseat(DAgent) = 0;
                                               seat2_interference =  seat2_interference+1;
                                            %If Column F and D are occupied
                                            elseif inseat(FAgent) ~= 0 && inseat(DAgent) ~= 0 && inseat(EAgent) == 0 && Coordinates(2,1,EAgent,t) > rowy(i) - pitchY/2 && Coordinates(2,1,EAgent,t) < rowy(i) + pitchY/2
                                               Coordinates(:,1,DAgent,t) = [widthX+smallx;rowy(i)-abs(life(i)-life2(i))*pitchY];
                                               blocknmoving(DAgent) = 1;
                                               blocknmoving(FAgent) = 0;
                                               inseat(DAgent) = 0;
                                               seat1_interference =  seat1_interference+1;
                                            end
                                       end
                                end
         
                        elseif  Coordinates(2,1,i,t) > rowy(i) + pitchY/2

                            
                                if Coordinates(1,1,i,t) < widthX
                                    Coordinates(1,1,i,t) = widthX + smallx;
                                elseif Coordinates(1,1,i,t) >  widthX + width_walk;
                                    Coordinates(1,1,i,t) = widthX + width_walk - smallx;
                                end
                                
                                Rie = [widthX+width_walk/2;rowy(i)] - Coordinates(:,:,i,t);
                                Desireddirections(:,:,i,t) = Rie./norm(Rie);
                                Forces(:,:,i,t) = Forces(:,:,i,t) + (1/T_i).*(Desiredspeed(i,t).*Desireddirections(:,:,i,t) - Velocities(:,:,i,t));
                                
                                r = diameter(i);
                                LiW =  Coordinates(:,:,i,t) - [widthX;Coordinates(2,:,i,t)];
                                RiW =  Coordinates(:,:,i,t) - [widthX+width_walk;Coordinates(2,:,i,t)];
                                dL = norm(LiW);
                                dR = norm(RiW);
                                
                                nLiW = LiW./dL;
                                nLiWt = [-nLiW(2);nLiW(1)];
                                nRiW = RiW./dR;
                                nRiWt = [-nRiW(2);nRiW(1)];
                    
                                if dR <= r
                                    Forces(:,:,i,t) = Forces(:,:,i,t) + RiW.*(U*exp(-dR/R)) + k*(r-dR).*nRiW - (K*(r-dR)*dot(Velocities(:,:,i,t),nRiWt)).*nRiWt;
                                else
                                    Forces(:,:,i,t) = Forces(:,:,i,t) + RiW.*(U*exp(-dR/R));
                                end
                                
                                if dL <= r
                                    Forces(:,:,i,t) = Forces(:,:,i,t) + LiW.*(U*exp(-dL/R)) + k*(r-dL).*nLiW - (K*(r-dL)*dot(Velocities(:,:,i,t),nLiWt)).*nLiWt;
                                else
                                    Forces(:,:,i,t) = Forces(:,:,i,t) + LiW.*(U*exp(-dL/R));
                                end         
                        end
                    end
    end
    for i = 1:crowdsize
            Accelerations(:,:,i,t) = Forces(:,:,i,t)./mass(i,1);
            Velocities(:,:,i,t) = Velocities(:,:,i,t) + 0.5.*Accelerations(:,:,i,t);
            if norm(Velocities(:,:,i,t)) > Maxspeed(i,t)
                arg = atan(Velocities(2,1,i,t)/Velocities(1,1,i,t));
                if Velocities(1,1,i,t) < 0
                    Velocities(:,:,i,t) = [-Maxspeed(i,t)*cos(arg);Maxspeed(i,t)*sin(arg)];
                else
                    Velocities(:,:,i,t) = [Maxspeed(i,t)*cos(arg);Maxspeed(i,t)*sin(arg)];
                end
            end
    end
    for i = 1:crowdsize
        if Coordinates(2,1,i,t) < pitch_total && Coordinates(1,1,i,t) >= widthX && Coordinates(1,1,i,t) <=  widthX + width_walk && release_time(i) < t
            if  inplane(i) == 0
               inplane(i) = t;
               count = count+1;
            end
        end
    end
    %Queing Management
    for i = 1:crowdsize
       if t < release_time(i)
         Coordinates(:,:,i,t) = [-3;pitch_total+2];
       end
    end
end


% Validate data
 for i = 1:crowdsize
     if inplane(i) > inseat(i)
         fail(i) = i;
     end
     if start_lug(i) > stop_lug(i)
         fail2(i) = i;
     end
 end 



% Calculation part
inplane(inplane <= 0) = inf;
inseat(inseat <= 0) = inf;

RTcal = zeros(size(release_time));
RT_change = size(release_time);
for i = 1:length(release_time)
RTcal(i) = sum(release_time==release_time(i));
RT_change = RT_change + RTcal(i);
end




% Time when the first passenger is on an airplane
start_t = min(inplane);
% Time when the first passenger is sitting
sit_t = max(inseat);
% Time for boarding
    if min(inseat) == 0;
        time = NaN;
    else
        time = sit_t - start_t;
    end

% Time for passenger to take a seat. No one wants to wait in queue!
for i = 1:crowdsize
    timetoseat(i) = inseat(i) - inplane(i);
end
    

% Simplify x and y coordinate notations
for t = 1:length(timelimit)
    for i = 1:crowdsize
        X(i,t) = Coordinates(1,1,i,t);
        Y(i,t) = Coordinates(2,1,i,t);
    end
end



Video = VideoWriter('final-random.avi');
open(Video)


%--- Plotting Section ---%
for t = 1:length(timelimit)
    hold off
    plot(X(:,t),Y(:,t),'g.','MarkerSize',20)
    axis ([-0.5 2*widthX+width_walk+0.5 -0.25 pitch_total+0.25])
    pbaspect([2*widthX+width_walk+1 pitch_total+5 1])
    %axis ([-5 3.8 0 23])
    %pbaspect([8.8 23 1])
    hold on
    plot(planrx,plany,'k','linewidth',2);
    plot(planx,plany,'k','linewidth',2);
    plot(planTopx,planTopy,'r','linewidth',2);
    plot(planBtmx,planBtmy,'k','linewidth',2);
    plot(planSpawnx,planSpawny,'r','linewidth',2);
    axis off;
    drawnow 
    %pause(0.01)
   
    text(-3,5,num2str(t))
    frame = getframe;
    writeVideo(Video,frame)
    
end

frame = getframe;
writeVideo(Video,frame);
close(Video);