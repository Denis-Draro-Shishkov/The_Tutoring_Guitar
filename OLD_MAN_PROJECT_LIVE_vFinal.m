% Denis Shishkov
% Jose Puerta
% Tutoring Guitar - Guitar Audio Processing 
% 05/05/17
%{
%% Recording Data
Fs = 24000;
pause(2)
TimeLength = 8; % Time taken in seconds
numbits = 24;
Toby_Demonstration_05_04 = audiorecorder(Fs,numbits,1,-1);
%pause(5)
disp('Start Playing.')
recordblocking(Toby_Demonstration_05_04,TimeLength);
disp('End of Recording')
Testing = getaudiodata(Toby_Demonstration_05_04);
Audio_File = Testing;
%}
%% If using saved recObj
Demonstration = Toby_Demonstration_05_04;%;%Camillo_Playing;%One_Note_147_PiezoElectric; %Six_Chord_82_110_147_196_247_329_22050Fs; %One_Note_147_PiezoElectric;
Fs = Demonstration.SampleRate;  
numbits = Demonstration.BitsPerSample;
Audio_File = getaudiodata(Demonstration);
TimeLength = length(Audio_File)/Fs;
%% Processing Data
%[y,Fs] = audioread('C:\Users\Draro\Documents\MATLAB\AudioTesting.wav');
%Audio_File = Please_Let_This_Work_PiezoElectric; %
%Fs = 44100;
%TimeLength = 8; % Time taken in seconds
tic
Freqs = linspace(-Fs/2,Fs/2,Fs);
NV = 48; % Higher Number voices per Octave gives a better frequency resolution at the cost of processing power. Can go up to 48 with the cwt function.
NO = 9; % Higher Number of Octaves increases the range of frequencies. NO of 8 gives a frequency range going from 79 Hz to 20000 Hz assuming a Fs of 44100 Hz. Probably should make this variable later based on other input criteria. Don't need samples above 4 kHz as all guitar harmonics are pretty much non-existent past that point.
%BeggingyouWork = upsample(Audio_File,2);
[wt2, f2] = cwt(real(Audio_File()),'bump',Fs,'VoicesPerOctave',NV,'NumOctaves',NO); % Applying cwt function
%[wt2, f2] = cwt(Audio_File(),'morse',Fs,'VoicesPerOctave',NV,'NumOctaves',NO,'TimeBandwidth',120); % Applying cwt function
FrequencyLimit_Low = length(f2);    
FrequencyLimit_High = 160; % Index of maximum applicable frequency (2,000 Hz)

Processed_Signal = wt2; %(FrequencyLimit_High:FrequencyLimit_Low,:); % Limiting the Frequency Band
Signal_Frequency = f2; %(FrequencyLimit_High:FrequencyLimit_Low);
Signal_Time = (1:length(wt2))/size(wt2,2)*TimeLength;
toc

%% Okay so real chord processing
%% Screwing around with frequencies
ActualGuitarNotes = [82 87 92 98 104 110 117 123 131 139 147 156 165 175 185 196 208 220 233 247 262 277 294 311 329 349 370 392 415 440 466 494 523 554 587 622 659 698 740 784 831 880 932 988 1047];
Guitar_Notes_Repmat = repmat(ActualGuitarNotes,length(Signal_Frequency),1);
Actual_Frequencies_Repmat = repmat(Signal_Frequency,1,length(ActualGuitarNotes));
[~, Actual_Index_Played] = min(abs(Actual_Frequencies_Repmat-Guitar_Notes_Repmat),[],1);
Used_Frequency = Signal_Frequency(Actual_Index_Played);
Signal_Of_Interest = Processed_Signal(Actual_Index_Played,:);
%% Trying to subtract a bit after the first detected bump so that we avoid the bullshit of multiple detected peaks blanking out the intermediate frequencies
Signal_Frequency = Used_Frequency;
Processed_Signal = Signal_Of_Interest;

%plot(abs(Signal_Of_Interest(26,:)),'b')
lag = 1000;
delay = lag/2/Fs; % Approximate delay from the smoothing
output = abs(Signal_Of_Interest);
for ii = (1:4)
    output = smoothts(output,'b',lag);
end
%% Subtract each frequency by a scaled version of those around it? Ughhhhhhhhhhhhhhhhhhhh UGLY
if ishandle(42)
close(42)
end
if ishandle(4242)
close(4242)
end

Output_Adjust = output; % Adjusting output for microphone attenuation at lower frequencies
Output_Adjust(1,:) = 0.19/0.067*1.5*output(1,:);
Output_Adjust(2,:) = 0.19/0.107*1.5*output(2,:);
Output_Adjust(3,:) = 0.0856/0.0297*1.5*output(3,:);
Output_Adjust(4,:) = 0.0725/0.0373*1.5*output(4,:);
Output_Adjust(5,:) = 1.25*1.5*output(5,:);% Used to be: Output_Adjust(5,:) = 0.0751/0.0482*1.5*output(5,:);
Output_Adjust(6,:) = 1/1*1.5*output(6,:);
Output_Adjust(7,:) = 1/1*1.25*output(7,:);
Output_Adjust(8,:) = 1/1*1.15*output(8,:);
Output_Adjust(9,:) = 1/1*1.05*output(9,:);
%{
figure(42)
output_2 = abs(Signal_Of_Interest);
for ii = 1:size(output_2,1)
    plot(linspace(0,TimeLength,length(Processed_Signal)),output_2(ii,:),'b')
    plot(linspace(0,TimeLength,length(Processed_Signal)),output_2(1,:),'r')
    [max_y, max_x] = max(output_2(ii,:));
    max_x = max_x/Fs;
    text(max_x,max_y,num2str(ii))   
    hold on;
end
title('Original Frequency Slices')
xlabel('Time (sec)')
ylabel('Frequency Intensity')
set(gca,'fontsize',18)

figure(4242)
for ii = 1:size(Output_Adjust,1)
    plot(linspace(0,TimeLength,length(Processed_Signal)),Output_Adjust(ii,:),'r')
    plot(linspace(0,TimeLength,length(Processed_Signal)),Output_Adjust(1,:),'b')
    [max_y, max_x] = max(Output_Adjust(ii,:));
    max_x = max_x/Fs;
    text(max_x,max_y,num2str(ii))
    hold on;
end
%}
title('Frequency Slices Adjusted -> Smoothed and Scaled')
xlabel('Time (sec)')
ylabel('Frequency Intensity')
set(gca,'fontsize',18)
%% Subtraction of Harmonics
output_new = zeros(length(ActualGuitarNotes)+1,size(wt2,2));
output_new(1,:) = Output_Adjust(1,:)*.8;
output_new(2:46,:) = Output_Adjust;

Output_HarmonicAdjust = zeros((size(Output_Adjust,1)+13),size(Output_Adjust,2));
Output_HarmonicAdjust(1,:) = Output_Adjust(1,:)*.8; % Generating one extra frequency to account for spectrum of the frequencies overlapping
Output_HarmonicAdjust(2:(size(Output_Adjust,1)+1),:) = Output_Adjust;
for ii = 13:size(Output_Adjust,1) 
    Output_HarmonicAdjust((ii),:) = output_new(ii,:) - output_new(ii-12,:)*1.7; % Subtraction of 12 because of octaves (12) 
    Output_HarmonicAdjust((ii),:) = smoothts(Output_HarmonicAdjust((ii),:),'b',lag);
end

for ii = 20:size(Output_Adjust,1) 
    Output_HarmonicAdjust((ii),:) = Output_HarmonicAdjust(ii,:) - output_new(ii-19,:)*1.4; % Subtraction of 19 because of third harmonic 
    Output_HarmonicAdjust((ii),:) = smoothts(Output_HarmonicAdjust((ii),:),'b',lag);
end

Output_HarmonicAdjust = Output_HarmonicAdjust(2:(size(Output_Adjust,1)+1),:);
%%
if ishandle(29)
close(29)
end
%{
figure(29)
for ii = 1:size(Output_HarmonicAdjust)
    %plot(Output_HarmonicAdjust(ii,:),'b')
    plot(linspace(0,TimeLength,length(Processed_Signal)),Output_HarmonicAdjust(ii,:),'b')   
    [max_y, max_x] = max(Output_HarmonicAdjust(ii,:));
    max_x = max_x/Fs;
    text(max_x,max_y,num2str(ii))
    hold on;
end
%}
title('Frequency Slices Adjusted -> Subtracted Harmonics and Smoothed')
xlabel('Time (sec)')
ylabel('Frequency Intensity')
set(gca,'fontsize',18)
%%
Pad_Zeros = zeros(1,size(Output_HarmonicAdjust,1));
Derivative1 = [Pad_Zeros' diff(abs(Output_HarmonicAdjust),1,2)]; 
Derivative2 = [Pad_Zeros', Pad_Zeros', diff(abs(Output_HarmonicAdjust),2,2)];

GetRidOf = zeros(size(Output_HarmonicAdjust,1),size(Output_HarmonicAdjust,2));
for aa = 1:size(Output_HarmonicAdjust)
    GetRidOf(aa,:) = ((Derivative1(aa,:) >= 5e-6)); %& (Derivative2(aa,:) > 5e-6));
end
%{
if ishandle(1)
close(1)
end
figure(1)
plot(Derivative1(3,:))
hold on;
plot(Derivative1(4,:))
%}
%%
if ishandle(966)
close(966)
end
%Changes = diff(abs(Processed_Signal),1,2); % Changes contains the difference between each index - i.e. rate of change
%Changes = abs(Processed_Signal);
%Changes = output;
Changes = Output_HarmonicAdjust.*GetRidOf;
%Changes = Output_HarmonicAdjust;
%LocationofAllPeaks = zeros(length(Audio_File),length(Signal_Frequency));
%IndsofAllPeaks = zeros(length(Audio_File),length(Signal_Frequency));
Freq_Cutoff = 0.01;
Global_Peaks_Ref = zeros(length(Signal_Frequency),length(Audio_File));

for ia = 1:length(Signal_Frequency)
    [Peaks_ia, Peaks_ia_Inds] = findpeaks(Changes(ia,:)); % Finds any peak rate of change for each frequency
    Global_Peaks_Ref(ia,Peaks_ia_Inds) = Peaks_ia; % Contains an array of peak intensities for each index. Each row is the indicies of the frequency.
end
Global_Peaks_Ref = Changes;
%{
figure(966)
for id = 1:length(Signal_Frequency)
    plot(linspace(0,TimeLength,length(Processed_Signal)),Global_Peaks_Ref(id,:),'b')   
    hold on;
end
%}
[Intensity_Values, Max_Peak_Time_Ind] = max(Global_Peaks_Ref,[],1);

Removed_Indicies = (Intensity_Values < Freq_Cutoff); %&& (Global_Peaks_Ref > 0)); Determines the indicies that we will be ignoring because of the too low peak.
Max_Peak_Time_Ind(Removed_Indicies) = 0; % Sets the removed indicies to 0. If index is 0, no LED is lit. Otherwise the LED corresponding to that specific frequency is lit. 


%% Afterwards Processing
%Notes = %Signal_Frequency(FreqIndex);
tic
Num_InDistinct_Notes = length(Max_Peak_Time_Ind); %length(Notes);3
ActualGuitarNotes = [82 87 92 98 104 110 117 123 131 139 147 156 165 175 185 196 208 220 233 247 262 277 294 311 329 349 370 392 415 440 466 494 523 554 587 622 659 698 740 784 831 880 932 988 1047 0];
FlippedGuitarNotes = fliplr(ActualGuitarNotes);

Notes_Compare = repmat(Max_Peak_Time_Ind,length(ActualGuitarNotes),1);
ActualGuitarNotes_Compare = repmat(ActualGuitarNotes,Num_InDistinct_Notes,1)';
[~, Actual_Index_Played] = min(abs(ActualGuitarNotes_Compare - Notes_Compare),[],1);

WORK = Max_Peak_Time_Ind; %ActualGuitarNotes(Actual_Index_Played);
toc
%% Loop that zeros out the false values that are too close to the actual played note -> This for loop isn't optimized. Takes forever to run.

[~,Time_Index] = find(WORK);
Output_Time_Index = Time_Index;
Index_Range = 150; % Value around which there shouldn't be any points - Index_Range/length(Audio_File)*TimeLength*1000 = Time Range getting zeroed
WORK_Range_Check = zeros(2*Index_Range+1,1);

for ii = 1:length(Time_Index)
    Compare_Time_Value = Time_Index(ii);
    Time_Index_Range = (Compare_Time_Value-Index_Range):(Compare_Time_Value+Index_Range); % Main Time Index +/- Index Range
    WORK_Range_Check = WORK(Time_Index_Range); % Range of values surrounding the identified one
    Valid_Indicies = find(WORK_Range_Check); 
    if (length(Valid_Indicies) > 1)
        Valid_Intensities = Intensity_Values(Time_Index_Range(Valid_Indicies)); % Add if statement to check for the notes being close together?
    
        [Max_Value, Max_Index] = max(Valid_Intensities);
    
        Removed_Indicies = Valid_Indicies;
        Removed_Indicies(Max_Index) = [];
        for jj = 1:length(Removed_Indicies)
            Output_Time_Index = Output_Time_Index(Output_Time_Index~=Compare_Time_Value+(Removed_Indicies(jj)-Index_Range-1)); % Deleting the Identified Index
        end
    end
end

% Output_Time_Index in the end has only the identified notes
%% Generating File that'll be imported into the Arduino
Frequency_Index = WORK(Output_Time_Index);
Time_Delay = diff(Output_Time_Index).*1000./Fs; % Saving the time delay in the format of ms (This is how the arduino is expecting the input)

Output_File_To_Arduino = zeros(length(Output_Time_Index),3);
Full_Time_Start = zeros(length(Output_Time_Index),1);
%Output_File_To_Arduino(1,3) = Output_Time_Index(1).*1000./Fs;
%Output_File_To----------------+_Arduino(2:end,3)= Time_Delay';
Full_Time_Start(1) = Output_Time_Index(1).*1000./Fs;
Full_Time_Start(2:end) = Time_Delay';

LED_Index = 1;
String_Num = 11; 
LoopAdjust = 0;
for ii = 1:length(Output_Time_Index) % Mother of all if statements. I actually can't think of a better way of doing this.... besides doing this 120 times. What the hell. 
    for jj = 1 % This loop only exists so that this beautiful code can be minimized in MATLAB
        if (Frequency_Index(ii) == 1) % 82
            LED_Index = 32; % Open Fret: Needs its own LED
            String_Num = 1;
        end
        if (Frequency_Index(ii) == 2)
            LED_Index = 1;
            String_Num = 1;
        end
        if (Frequency_Index(ii) == 3)
            LED_Index = 3;
            String_Num = 1;
        end 
        if (Frequency_Index(ii) == 4)
            LED_Index = 7;
            String_Num = 1;
        end
        if (Frequency_Index(ii) == 5)
            LED_Index = 5;
            String_Num = 1;
        end
        if (Frequency_Index(ii) == 6)
            LED_Index = 32;  % Open Fret: Needs its own LED - Problem with LEDs, so the code parameters are different
            String_Num = 2;
        end      
        if (Frequency_Index(ii) == 7)
            LED_Index = 1;
            String_Num = 2;
        end
        if (Frequency_Index(ii) == 8)
            LED_Index = 3;
            String_Num = 2;
        end
        if (Frequency_Index(ii) == 9)
            LED_Index = 7;
            String_Num = 2;
        end 
        if (Frequency_Index(ii) == 10)
            LED_Index = 5;
            String_Num = 2;
        end
        if (Frequency_Index(ii) == 11)
            LED_Index = 32;  % Open Fret: Needs its own LED
            String_Num = 3;
        end
        if (Frequency_Index(ii) == 12)
            LED_Index = 1;
            String_Num = 3;
        end
        if (Frequency_Index(ii) == 13)
            LED_Index = 3;
            String_Num = 3;
        end
        if (Frequency_Index(ii) == 14)
            LED_Index = 7;
            String_Num = 3;
        end
        if (Frequency_Index(ii) == 15)
            LED_Index = 5;
            String_Num = 3;
        end 
        if (Frequency_Index(ii) == 16)
            LED_Index = 32; % Open Fret: Needs its own LED
            String_Num = 4;
        end
        if (Frequency_Index(ii) == 17)
            LED_Index = 1;
            String_Num = 4;
        end
        if (Frequency_Index(ii) == 18)
            LED_Index = 3;
            String_Num = 4;
        end      
        if (Frequency_Index(ii) == 19)
            LED_Index = 7;
            String_Num = 4;
        end
        if (Frequency_Index(ii) == 20)
            LED_Index = 32; % Open Fret: Needs its own LED
            String_Num = 5;
        end
        if (Frequency_Index(ii) == 21)
            LED_Index = 1;
            String_Num = 5;
        end 
        if (Frequency_Index(ii) == 22)
            LED_Index = 3;
            String_Num = 5;
        end
        if (Frequency_Index(ii) == 23)
            LED_Index = 7;
            String_Num = 5;
        end
        if (Frequency_Index(ii) == 24)
            LED_Index = 5;
            String_Num = 5;
        end    
        if (Frequency_Index(ii) == 25)
            LED_Index = 32;  % Open Fret: Needs its own LED
            String_Num = 6;
        end 
        if (Frequency_Index(ii) == 26)
            LED_Index = 1;
            String_Num = 6;
        end
        if (Frequency_Index(ii) == 27)
            LED_Index = 3;
            String_Num = 6;
        end
        if (Frequency_Index(ii) == 28)
            LED_Index = 7;
            String_Num = 6;
        end      
        if (Frequency_Index(ii) == 29)
            LED_Index = 5;
            String_Num = 6;
        end
        if (Frequency_Index(ii) == 30)
            LED_Index = 8;
            String_Num = 6;
        end
        if (Frequency_Index(ii) == 31)
            LED_Index = 16;
            String_Num = 6;
        end 
    end % Ze Mother of all if statements
    Output_File_To_Arduino(ii+LoopAdjust,1) = LED_Index; %ActualGuitarNotes(Frequency_Index(ii));%LED_Index;
    Output_File_To_Arduino(ii+LoopAdjust,2) = (String_Num-1); % Minus 1 because indexing in arduino IS DUMB
    Output_File_To_Arduino(ii+LoopAdjust,3) = Full_Time_Start(ii);
    % Inserting Function that would make it possible to flash the LEDs
    % after they are turned on.    
end

format long g
Output_File_To_Arduino % Before LED flash Code

%% Inserting Flash LED code - WORKS
Full_Time_Final = Full_Time_Start';

Final_LED_Index = Output_File_To_Arduino(:,1)';
Final_String_Num = Output_File_To_Arduino(:,2)';
Final_Time_Index = Output_File_To_Arduino(:,3)';
Full_Time_Final = Final_Time_Index;

LED_ON_TIME = 500; % Blink LEDs after 500ms. 
Insertion = 0; % Counts how many elements have been inserted. 
for ii = 1:(length(Output_Time_Index)-1)
    Index = length(Output_Time_Index) - ii; % Going backwards through the loop
    String_Num = Output_File_To_Arduino(Index,2);
    LED_Index = Output_File_To_Arduino(Index,1);
    
    %End_Index = find(Full_Time_Final(Index:end) < LED_ON_TIME)
    Temp_Time = LED_ON_TIME;
    if (Index == length(Output_Time_Index)-1) % Final Element %isempty(End_Index) % Don't need to adjust the indicies anymore. 
        %test = ii
        Final_String_Num(Index+1) = String_Num; % Minus 1 because indexing in arduino IS DUMB
        Final_Time_Index(Index+1) = Temp_Time;
        Final_Time_Final(Index+1) = Temp_Time;
        Final_LED_Index(Index+1) = LED_Index + 200;
        
        % Delete 500 from next peak.
    else
        End_Index = length(Final_Time_Index);
        
        for jj = (Index):(End_Index-1)
            if ((Temp_Time - Final_Time_Index(jj+1)) < 0)
                break; 
            else
                Temp_Time = Temp_Time - Final_Time_Index(jj+1); % This + 1 seems to make all the difference. NEED TO CHECK!
            end
        end
        
        % took away + Insertion
        Final_String_Num = [Final_String_Num(1:jj) String_Num Final_String_Num(jj+1:end)]; % Insertion to say when to flash
        Final_Time_Index(jj+1) = Final_Time_Index(jj+1) - Temp_Time;
        Final_Time_Index = [Final_Time_Index(1:jj) Temp_Time Final_Time_Index(jj+1:end)];
        Final_LED_Index = [Final_LED_Index(1:jj) (LED_Index + 200) Final_LED_Index(jj+1:end)];
        Full_Time_Final = [Full_Time_Final(1:jj) Temp_Time Full_Time_Final(jj+1:end)];
    end
    
    Insertion = Insertion + 1;
    
end
Output_File_To_Arduino_Final = zeros(length(Final_LED_Index),3);

Final_Time_Index(1:(end-1)) = Final_Time_Index(2:(end));

Output_File_To_Arduino_Final(:,1) = Final_LED_Index';
Output_File_To_Arduino_Final(:,2) = Final_String_Num';
Output_File_To_Arduino_Final(:,3) = Final_Time_Index';
format long g
%Output_File_To_Arduino 
Output_File_To_Arduino_Final

%%
Output_File_To_Arduino_Adjusted = Output_File_To_Arduino_Final';
fileID = fopen('SeniorProjectTextFiles/Code_To_Guitar.txt','w');
fprintf(fileID,'  {%.0f,%.0f,%.0f},\n',Output_File_To_Arduino_Adjusted); % To get higher accuracy, you can split up the arduino code into millidelay and microdelay. Is it worth it?
fclose(fileID);
