clearvars;
close all;

% This code will run a simulation (or multiple depending on what is
% determined by "Iterations") of proteins binding and unbinding from a DNA
% lattice. When a protein binds to a given location a timer will begin and
% continue timing until the protein unbinds. Once the protein unbinds, the
% timer stops and the total dwell time is recorded. This is performed for
% each protein that binds AND unbinds from the lattice (proteins left on at
% the end will not be counted as no sufficient timing data can be
% recorded). The dwell times will then ue placed in a histogram to compare
% to KERA data.

N = 8660;   %length of DNA lattice
n = 3;  %length of a monomer
w = 1;  %cooperativity parameter
L_Total = 2;    %total concentration of RAD51
k_on = 1;   %kinetic rate constant
k_off = 1;
Ratio = 1;   %Percentage of solution which is monomers (0 to 1)
Iterations = 5;    %number of iterations of the simulation (number of growth profiles generated)

UncoveredLength = 0.34; %length of a DNA nt without RAD51 bound to it (according to van der Heijden paper) - nm
CoveredLength = 0.51;   %length of a DNA nt where RAD51 is bound - nm

minEvents = 1000;   %minimum number of events simulated with each iteration (used for ending at equilibrium)

%Memory Allocation
EventFractions = zeros(Iterations,7);
FracCover = zeros(Iterations,minEvents);
DNA_Length = zeros(Iterations,minEvents);
t = zeros(Iterations,minEvents);
Max_Time = zeros(1,Iterations);
Equilibrium_Coverage = zeros(1,Iterations);
Dwell_Times = [];   %empty array that will be base of dwell times (times will be added here when a protein unbinds)

L_Monomer = Ratio*L_Total;        %Concentration of monomer RAD51
L_Dimer = (1-Ratio)*L_Total;    %Concentration of dimer RAD51

for Loops = 1:Iterations
    DNA = zeros(1,N);
    BoundAtSpot = zeros(1,N);   %records where monomers are bound on lattice
    Protein_Timers = zeros(1,N); %timer at each location for timing how long proteins remain on the lattice
       
    %Memroy Allocations
    Populations = zeros(minEvents,7);
    a = zeros(minEvents,7);
    Probabilities = zeros(minEvents,7);
    FiringAmounts = zeros(minEvents,7);
    dt = zeros(1,minEvents);
    j = zeros(1,minEvents);
    Location_History = zeros(7,minEvents);

    Equilibrium = 0;
    Events = 0;
    while max(t(Loops,:)) < 1.25    %time value the system runs to
        Events = Events+1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        searcH_vars = {'Left','Right','Gap_Size','Left_Available_M','Right_Available_M','Left_Available_D','Right_Available_D','Right_BindingSite_M','Right_BindingSite_D','Gaps_L2_M','Gaps_R2_M','Gaps_L2_D','Gaps_R2_D','Gap_Size2_M','Gap_Size2_D','Gap_SizeI_M','Gap_SizeI_D','Left_I_M','Left_I_D','Isolated_M','SinglyContiguous_M','Doubly_Contiguous_M','Isolated_D','SinglyContiguous_D','DoublyContiguous_D'};
        clear searcH_vars;

        Left = find(diff([1 DNA 1]) == -1);
        Right = find(diff([1 DNA 1]) == 1)-1;
        Gap_Size = Right-Left+1;

        Left_Available_M = Left(Gap_Size >= n);
        Right_Available_M = Right(Gap_Size >= n);
        Left_Available_D = Left(Gap_Size >= 2*n);
        Right_Available_D = Right(Gap_Size >= 2*n);
        Right_BindingSite_M = Right_Available_M-(n-1);
        Right_BindingSite_D = Right_Available_D-(2*n-1);

        %Doubly Contiguous Searches
        Doubly_Contiguous_M = Left_Available_M(Left_Available_M == Right_BindingSite_M & 1 < Left_Available_M & Left_Available_M < N-(n-1));
        Doubly_Contiguous_D = Left_Available_D(Left_Available_D == Right_BindingSite_D & 1 < Left_Available_D & Left_Available_D < N-(2*n-1));

        %Singly Contiguous Searches
        Singly_Contiguous_M = unique([Left_Available_M(~ismember(Left_Available_M,Doubly_Contiguous_M)),Right_Available_M(~ismember(Right_Available_M-(n-1),Doubly_Contiguous_M))-(n-1)]);
        Singly_Contiguous_M(Singly_Contiguous_M == N-(n-1) | Singly_Contiguous_M == 1) = [];
        Singly_Contiguous_D = unique([Left_Available_D(~ismember(Left_Available_D,Doubly_Contiguous_D)),Right_Available_D(~ismember(Right_Available_D-(2*n-1),Doubly_Contiguous_D))-(2*n-1)]);
        Singly_Contiguous_D(Singly_Contiguous_D == N-(2*n-1) | Singly_Contiguous_D == 1) = [];

        %Isolated Searches
        Gaps_L2_M = Left(~ismember(Left,Doubly_Contiguous_M));
        Gaps_L2_D = Left(~ismember(Left,Doubly_Contiguous_M) & ~ismember(Left,Doubly_Contiguous_D));
        Gaps_R2_M = Right(~ismember(Right,Doubly_Contiguous_M+(n-1)));
        Gaps_R2_D = Right(~ismember(Right,Doubly_Contiguous_M+(n-1)) & ~ismember(Right,Doubly_Contiguous_D+(2*n-1)));
        Gap_Size2_M = Gaps_R2_M-Gaps_L2_M+1;
        Gap_Size2_D = Gaps_R2_D-Gaps_L2_D+1;

        Gap_SizeI_M = Gap_Size2_M(Gap_Size2_M > n); %size of gaps where isolated binding can happen
        Gap_SizeI_D = Gap_Size2_D(Gap_Size2_D > 2*n);
        Left_I_M = Gaps_L2_M(ismember(Gap_Size2_M,Gap_SizeI_M));    %left position of gaps where isolated binding can happen
        Left_I_D = Gaps_L2_D(ismember(Gap_Size2_D,Gap_SizeI_D));

        Isolated_M = zeros(1,sum(Gap_SizeI_M)-((n+1)*numel(Gap_SizeI_M))+logical(ismember(N,Left_I_M+Gap_SizeI_M-1))+logical(sum(DNA) == 0));   %memory allocation
        Isolated_D = zeros(1,sum(Gap_SizeI_D)-((2*n+1)*numel(Gap_SizeI_D))+logical(ismember(N,Left_I_D+Gap_SizeI_D-1))+logical(sum(DNA) == 0));
        for i = 1:numel(Left_I_M)
            Isolated_M(find(Isolated_M == 0,1):find(Isolated_M == 0,1)+Gap_SizeI_M(i)-(n+1)-1+logical(N == Left_I_M(i)+Gap_SizeI_M(i)-1)+logical(Left_I_M(i) == 1)) = Left_I_M(i)+(1-logical(Left_I_M(i) == 1):Gap_SizeI_M(i)-(n+1)+logical(N == Left_I_M(i)+Gap_SizeI_M(i)-1));
        end
        for k = 1:numel(Left_I_D)
            Isolated_D(find(Isolated_D == 0,1):find(Isolated_D == 0,1)+Gap_SizeI_D(k)-(2*n+1)-1+logical(N == Left_I_D(k)+Gap_SizeI_D(k)-1)+logical(Left_I_D(k) == 1)) = Left_I_D(k)+(1-logical(Left_I_D(k) == 1):Gap_SizeI_D(k)-(2*n+1)+logical(N == Left_I_D(k)+Gap_SizeI_D(k)-1));
        end

        %Population Numbers
        xB_IM = length(Isolated_M);
        xB_SCM = length(Singly_Contiguous_M);
        xB_DCM = length(Doubly_Contiguous_M);
        xB_ID = length(Isolated_D);
        xB_SCD = length(Singly_Contiguous_D);
        xB_DCD = length(Doubly_Contiguous_D);
        xAB = sum(DNA)/n;

        Populations(Events,:) = [xB_IM,xB_SCM,xB_DCM,xB_ID,xB_SCD,xB_DCD,xAB];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        a(Events,:) = Populations(Events,:).*[k_on*L_Monomer,k_on*L_Monomer*w,k_on*L_Monomer*(w^2),k_on*L_Dimer,k_on*L_Dimer*w,k_on*L_Dimer*(w^2),k_off];   %calculation of all propensity functions
        Probabilities(Events,:) = a(Events,:)./sum(a(Events,:));

        r = [rand,rand,rand,rand,rand,rand,rand];
        tau = (1./a(Events,:)).*log(1./r);  %time step to next firing of each reaction
        FiringAmounts(Events,:) = a(Events,:).*tau; %average amount of firing of each reaction in the time step
        dt(Events) = min(tau);  %actual time step that occurs (minimum of tau)
        j(Events) = find(tau == min(tau));  %determination of which even occurs

        Protein_Timers(BoundAtSpot == 1) = Protein_Timers(BoundAtSpot == 1) + dt(Events);   %advances timer wherever a protein is bound
        
        if j(Events) == 1        %Isolated monomer binding
            Bind_Spot_I_M = Isolated_M(randi(xB_IM));
            DNA(Bind_Spot_I_M:Bind_Spot_I_M+(n-1)) = 1;
            Location_History(1,Events) = Bind_Spot_I_M;
            BoundAtSpot(Bind_Spot_I_M) = 1;
        elseif j(Events) == 2    %Singly-Contiguous monomer binding
            Bind_Spot_SC_M = Singly_Contiguous_M(randi(xB_SCM));
            DNA(Bind_Spot_SC_M:Bind_Spot_SC_M+(n-1)) = 1;
            Location_History(2,Events) = Bind_Spot_SC_M;
            BoundAtSpot(Bind_Spot_SC_M) = 1;
        elseif j(Events) == 3    %Doubly-Contiguous monomer binding
            Bind_Spot_DC_M = Doubly_Contiguous_M(randi(xB_DCM));
            DNA(Bind_Spot_DC_M:Bind_Spot_DC_M+(n-1)) = 1;
            Location_History(3,Events) = Bind_Spot_DC_M;
            BoundAtSpot(Bind_Spot_DC_M) = 1;
        elseif j(Events) == 4    %Isolated dimer binding
            Bind_Spot_I_D = Isolated_D(randi(xB_ID));
            DNA(Bind_Spot_I_D:Bind_Spot_I_D+(2*n-1)) = 1;
            Location_History(4,Events) = Bind_Spot_I_D;
            BoundAtSpot(Bind_Spot_I_D) = 1;
            BoundAtSpot(Bind_Spot_I_D+n) = 1;
        elseif j(Events) == 5    %Singly-Contiguous dimer binding
            Bind_Spot_SC_D = Singly_Contiguous_D(randi(xB_SCD));
            DNA(Bind_Spot_SC_D:Bind_Spot_SC_D+(2*n-1)) = 1;
            Location_History(5,Events) = Bind_Spot_SC_D;
            BoundAtSpot(Bind_Spot_SC_D) = 1;
            BoundAtSpot(Bind_Spot_SC_D+n) = 1;
        elseif j(Events) == 6    %Doubly-Contiguous dimer binding
            Bind_Spot_DC_D = Doubly_Contiguous_D(randi(xB_DCD));
            DNA(Bind_Spot_DC_D:Bind_Spot_DC_D+(2*n-1)) = 1;
            Location_History(6,Events) = Bind_Spot_DC_D;
            BoundAtSpot(Bind_Spot_DC_D) = 1;
            BoundAtSpot(Bind_Spot_DC_D+n) = 1;
        elseif j(Events) == 7    %Monomer unbinding
            Bound_Locations = find(BoundAtSpot == 1);
            Unbind_Spot = Bound_Locations(randi(length(Bound_Locations)));
            DNA(Unbind_Spot:Unbind_Spot+(n-1)) = 0;
            Location_History(7,Events) = Unbind_Spot;
            BoundAtSpot(Unbind_Spot) = 0;
            Dwell_Times = [Dwell_Times,Protein_Timers(Unbind_Spot)+dt(Events)]; %adds time on lattice to dwell time array; additional dt must be added since that's how long it took for an unbinding event to happen since the last event
            Protein_Timers(Unbind_Spot) = 0;    %resets Protein_Timers at the location of unbinding
        end
        
        FracCover(Loops,Events+1) = sum(DNA)/N; %calculates saturation of DNA lattice
        t(Loops,Events+1) = t(Loops,Events)+dt(Events); %updates time value
        
        DNA_Length(Loops,Events+1) = (numel(find(DNA==1))*CoveredLength)+(numel(find(DNA==0))*UncoveredLength); %calculated length of DNA strand based on lengths given previously
        
        if ~isempty(find(BoundAtSpot == 1 & DNA ~= 1, 1))   %error where location is recorded as having a bound protein but it isn't actually present
            disp(['STOP - ERROR AT ', num2str(find(BoundAtSpot == 1 & DNA ~= 1, 1))]);
            flag = 1;
            break
        end
    end

    EventFractions(Loops,:) = [numel(find(j==1)),numel(find(j==2)),numel(find(j==3)),numel(find(j==4)),numel(find(j==5)),numel(find(j==6)),numel(find(j==7))]./Events;  %fraction of how often each event occured
    Max_Time(Loops) = max(t(Loops,:));  %absolute maximum time across all iterations
   
    figure(1);  %plot of growth profiles (subplot 1 = saturation, subplot 2 = length {van der Heijden})
%     subplot(2,1,1);
    hold on;
    scatter(t(Loops,:),FracCover(Loops,:),1,'b','filled','HandleVisibility','off');
    ylabel('Fractional Coverage');
    xlim([0 1.25]);
    ylim([0 1]);
    title('Saturation of DNA Lattice');
    box on;
%     subplot(2,1,2);
%     hold on;
%     scatter(t(Loops,:),DNA_Length(Loops,:)/1000,1,Colors(kon_Values == w,:),'filled','HandleVisibility','off');
%     ylabel('Length (\mum)');
    xlabel('Time, t');
%     xlim([0 1.25]);
%     ylim([N*UncoveredLength/1000 N*CoveredLength/1000]);
%     title('Length of DNA Strand');
%     box on;    

    if flag == 1
        break
    end 
end

figure(2);  %histogram of dwell times of proteins on lattice
histfit(Dwell_Times,50,'exponential');
hold on;
xlabel('Dwell Times');
xlim([0 inf]);
ylabel('Count');
title('Dwell Time Histogram');
box on;