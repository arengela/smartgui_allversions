function varargout = smartGUI_STFT(varargin)
% smartGUI2 MATLAB code for smartGUI2.fig
%      smartGUI2, by itself, creates a new smartGUI2 or raises the existing
%      singleton*.
%
%      H = smartGUI2 returns the handle to a new smartGUI2 or the handle to
%      the existing singleton*.
%
%      smartGUI2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in smartGUI2.M with the given input arguments.
%
%      smartGUI2('Property','Value',...) creates a new smartGUI2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before smartGUI2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to smartGUI2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help smartGUI2

% Last Modified by GUIDE v2.5 15-May-2013 14:12:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
%uicontrol('Style','pushbutton','String','collect','Position',[10 10 20 20])
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @smartGUI_STFT_OpeningFcn, ...
    'gui_OutputFcn',  @smartGUI_STFT_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before smartGUI2 is made visible.
function smartGUI_STFT_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to smartGUI2 (see VARARGIN)

% Choose default command line output for smartGUI2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes smartGUI2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = smartGUI_STFT_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in start_flag.
function start_flag_Callback(hObject, eventdata, handles)
% hObject    handle to start_flag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Set variables
if isempty(get(handles.htkFilePath,'String'))
    simData=0;
else
    simData=1;
end

if simData==0
    sampling_rate = 400;
    number_of_channels = 64;
    ANnumber_of_channels = 4;
    number_of_sec2read = 2;
    dANproj_name =  {'Amp.dmata'};
    sANproj_name =  {'Amp.smata'};
    integration_time = .1;
    avgEvent_freqs2plot = 1:7;
    singleEvent_freqs2plot = 1:7;
    num_freq_bands = 1;
    window_around_event = 500;
    window_around_event_ms=1000;
    constant = []; posInNewData_AfterCAR = [];
    ANconstant = []; ANposInNewData = [];
else
    sampling_rate = 400;
    number_of_channels = 16;
    ANnumber_of_channels = 4;
    number_of_sec2read = 1;
    dANproj_name =  {'Amp.dmata'};
    sANproj_name =  {'Amp.smata'};
    integration_time = .1;
    avgEvent_freqs2plot = 1:17;
    singleEvent_freqs2plot = avgEvent_freqs2plot;
    num_freq_bands = length(avgEvent_freqs2plot);
    window_around_event = 50;
    window_around_event_ms = 1500;
     constant = [];posInNewData_AfterCAR=[];
      ANconstant = []; ANposInNewData = [];
end

% Flags
event_flag = 1;  % Average plot
single_event_flag = 1;
single_stacked_flag = 1;
ave_spec_flag=1;

% Retrieve from GUI
time2collectBaseline = str2num(get(handles.time2collectBaseline,'String'));
desired_freq_plot = max(avgEvent_freqs2plot) ; %str2num(get(handles.desired_freq_plot,'String'));
freq_band_singlestacked =  max(avgEvent_freqs2plot)-1; %str2num(get(handles.freq_band_singlestacked,'String'));
desired_ANchan = str2num(get(handles.desired_ANchan,'String'));
threshold = str2num(get(handles.threshold,'String'));
subj_id = get(handles.subj_id,'String');
elec_tmp = get(handles.get_total_electrodes,'String');
number_of_electrodes_total=str2num(elec_tmp{get(handles.get_total_electrodes,'Value')});
number_of_analog=length(desired_ANchan);


if  number_of_electrodes_total == 64%get(handles.get_total_electrodes,'Value')==1
    dproj_name = {'Amp.dmat1'};
    sproj_name = {'Amp.smat1'};
    to_plot_grid=8;
elseif number_of_electrodes_total == 256;
    dproj_name = {'Amp.dmat1';'Amp.dmat2';'Amp.dmat3';'Amp.dmat4'};
    sproj_name = {'Amp.smat1';'Amp.smat2';'Amp.smat3';'Amp.smat4'};
    to_plot_grid=16;
end

% Calculate
points_needed4baseline = sampling_rate*time2collectBaseline*number_of_channels *num_freq_bands;
number_of_points2read = sampling_rate*number_of_sec2read*number_of_channels*num_freq_bands;
ANnumber_of_points2read = sampling_rate*number_of_sec2read*ANnumber_of_channels;
num_avgEvent_freqs2plot = length(avgEvent_freqs2plot);
num_singleEvent_freqs2plot = length(singleEvent_freqs2plot);


%%
pth = pwd;
addpath([pth filesep 'movie']);
rehash path;


%%
bufferCounter = zeros(1,size(dproj_name,1));
ANbufferCounter = 0;
matlab_loopCounter = 0;
num_samples=0; %used to calculate faster average
last_used_ind=0;  %used to calculate faster average
eventRelatedAvg=[];


%buffer of 5 minutes
if event_flag ==1 || single_event_flag==1
    DataAfterCAR=zeros(number_of_electrodes_total, num_freq_bands, 4*60*sampling_rate);
    ANNewData_finalMAT=zeros(ANnumber_of_channels, 4*60*sampling_rate);
end
try
    DA = actxcontrol('TDevAcc.X');% Loads ActiveX methods
catch
    DA=[];
end
handles.DA = DA;

guidata(hObject, handles);

ii=1;

last_used_ind=-1;%%
good_event_count=0;
desired_ANchan = desired_ANchan(1);
try
    desired_ANchan2 = desired_ANchan(2);
catch
    desired_ANchan2=[];
end
% ADD TO GUI
%threshold = .5;
threshold2 = .7; % ADD TO GUI
num_events = 1;
indLastEvent = 0;
eventIndices = zeros(2,50); % 50 = max number of events expected
prev_num_events=0;
new_plot_flag=0;
finalPos=0;
finalPosA=0;
ecogFinalPos=0;
newEventFlag=0;
lastTime=0;
ecogHold.data=[];
ecogHold2.data=[];
specTimes2=[];
lastRawIdx=0;
nextRawIdx=1;
ecogAll.data=[];
ecogAll.sampFreq=sampling_rate;
if simData==1 || DA.ConnectServer('Local')>0  %Checks to see if there is an OpenProject and Workbench load
    DA.SetSysMode(3); pause(5)
    %    if DA.SetSysMode(2)==1%Starts OpenEx project running and acquiring data. 1 = Standby, 2 = Preview, !!MUST CHANGE TO 3 WHEN RECORDING
    if simData==1 || DA.GetSysMode>1 % Checks to see if the Project is running
        if simData==1 | not(isnan(DA.ReadTargetVEX(dproj_name{1},0,10,'F32','F64')))   %checks to see if it loaded the correct project
            try
                for i=1:length(dproj_name)
                    BufferSize(i)=DA.GetTargetSize(dproj_name{i});% returns the Buffer size (bookkeeping)
                end
                if event_flag ==1 || single_event_flag==1
                    ANBufferSize=DA.GetTargetSize(dANproj_name{1});% returns the Analog Buffer size (bookkeeping)
                    ANoldAcqPos=ANBufferSize;
                end
                oldAcqPos=BufferSize;
            end
            
            %if DA.GetSysMode=0, then means tdt is not recording or prj not working
            
            
            tic
            %%
            %Calculate baseline
            if ~isempty(get(handles.inputBaselineFile,'String')) & simData~=1
                load(get(handles.inputBaselineFile,'String'));
                medians = squeeze(medians(:,:,1));
            else
                if isempty(get(handles.inputBaselineFile,'String')) & simData~=1
                    [baseline,averages, stdevs, medians,logBaselineDataMat]=calculate_baseline(time2collectBaseline, ...
                        sproj_name,dproj_name, points_needed4baseline, BufferSize, DA, number_of_channels, ...
                        num_freq_bands, sampling_rate,number_of_sec2read);
                    pause(1)
                    uisave({'averages' 'stdevs' 'medians'}, ['baseline_stats_' subj_id])
                    medians = squeeze(medians(:,:,1));
                elseif simData==1
                    baseline=loadHTKBaseline(handles,time2collectBaseline);
                end
                %badChannels=detectBadChannels(baseline.data);
                badChannels=[];
                baseline=subtractCAR_16ChanBlocks(baseline,badChannels,baseline);
                
                baselineSTFT=calcSTFT(baseline.data,baseline.sampFreq,[]);
                sampPerSecond=size(baselineSTFT.data,3)/time2collectBaseline;
                window_around_event=round((window_around_event_ms/1000)*sampPerSecond);
                avgEvent_freqs2plot = 1:size(baselineSTFT.data,2);
                singleEvent_freqs2plot = avgEvent_freqs2plot;
                num_avgEvent_freqs2plot = length(avgEvent_freqs2plot);
                num_singleEvent_freqs2plot = length(singleEvent_freqs2plot);
                num_freq_bands = length(avgEvent_freqs2plot);
                allStacked=zeros(number_of_electrodes_total,50,num_freq_bands,window_around_event+1);
                desired_freq_plot=find(baselineSTFT.freq>70 & baselineSTFT.freq<150);
                
                medians=median(baselineSTFT.data,3);
                medians=repmat(medians,[1 1 window_around_event+1]);
                averages=mean(baselineSTFT.data,3);
                averages=repmat(averages,[1 1 window_around_event+1]);
                stdevs=std(baselineSTFT.data,[],3);
                stdevs=repmat(stdevs,[1 1 window_around_event+1]);
                
                 DataAfterCAR=zeros(number_of_electrodes_total, num_freq_bands, 4*60*sampling_rate);
                ANNewData_finalMAT=zeros(ANnumber_of_channels, 4*60*sampling_rate);
                
                handles=initiatePlots(handles,num_avgEvent_freqs2plot,num_singleEvent_freqs2plot,number_of_electrodes_total,...
                    window_around_event,to_plot_grid,desired_freq_plot,number_of_analog,sampPerSecond,allStacked)
            end
            
            
            
            profile on
            
            while (simData==1 || DA.GetSysMode>1) %& good_event_count<11
                %%
                %Retrieve data
                % find place in TDT buffer, check for buffer wrap
                % around, find place in MATLAB buffer
                if simData~=1;
                    [ecog,analog,LogNewData,ANNewData,index,logGdat,analog_dat,ANAcqPos,AcqPos,bufferCounter, ...
                        constant, ANconstant, posInNewData_AfterCAR, ANposInNewData, number_of_points2read, ANnumber_of_points2read]=getTDTData(handles,sampling_rate,ii,matlab_loopCounter,bufferCounter,ANbufferCounter,...
                        number_of_points2read, ANnumber_of_points2read,sproj_name,dproj_name,sANproj_name,dANproj_name,...
                        BufferSize,ANBufferSize, oldAcqPos, ANoldAcqPos,number_of_channels,ANnumber_of_channels, constant,ANconstant,...
                        posInNewData_AfterCAR, ANposInNewData);
                    %ANNewDataMAT=reshape(ANNewData,ANnumber_of_channels, ANnumber_of_points2read/ANnumber_of_channels);
                    ecogAll.data=[ecogAll.data ecog.data];
                else
                    %if ~exist('ecog') | newEventFlag==0
                    if matlab_loopCounter==0
                        [ecogAll1,analogAll]=loadHTKData(handles);
                    end
                    
                    [ecog,analog,lastTime]=grabDataBuffer(ecogAll1,analogAll,lastTime);                    

                    %end
                    
                    %LogNewData=ecogSTFT.data-repmat(median(baselineSTFT.data,3),[1 1 size(ecogSTFT.data,3)]);
                    %[LogNewData,ANNewData,index,logGdat,analog_dat,ANAcqPos,AcqPos,bufferCounter]=getSimData(handles,time2collectBaseline,sampling_rate,number_of_sec2read,ii,matlab_loopCounter,bufferCounter);
                end
                
                ecog=subtractCAR_16ChanBlocks(ecog,badChannels,ecog);
                ecogAll.data(:,ecogFinalPos+1:ecogFinalPos+size(ecog.data,2))=ecog.data;
                ecogFinalPos=ecogFinalPos+size(ecog.data,2);

                %ecog=applyLineNoiseNotch_60HzHarmonics(ecog,baseline.sampFreq);
                [ecogSTFT,nextRawIdx]=calcSTFT(ecogAll.data(:,nextRawIdx:ecogFinalPos),ecogAll.sampFreq,nextRawIdx);

                %ecogSTFT=ecog;
                %ecogSTFT.data=permute(ecogSTFT.data,[1,3,2]);
                ANNewData=analog.data;
               ANNewDataMAT=ANNewData;

                LogNewData=ecogSTFT.data;
                %%
                %Plot continual plot
                %continual plotting
                runningAverage = squeeze(mean(LogNewData(:, desired_freq_plot,:),3));
                
                %runningAverage = squeeze(mean(LogNewData(:, desired_freq_plot, (end-integration_time*sampling_rate):end),3));
                ZscoreNewData=mean((runningAverage-averages(:,desired_freq_plot,1))./stdevs(:,desired_freq_plot,1),2);
                to_plot=real(reshape(ZscoreNewData,to_plot_grid,to_plot_grid))';
                handles=plotContinual(handles,'update',[],[],to_plot);
                %%ter
                %Plot event coun
                % Place data into MATLAB buffers
                % finalPos = posInNewData_AfterCAR+sampling_rate*number_of_sec2read-1;
                %finalPos = posInNewData_AfterCAR + number_of_points2read/(num_freq_bands*number_of_channels) - 1;
                %ANfinalPos = ANposInNewData + ANnumber_of_points2read/ANnumber_of_channels - 1;
%                 ecogHold.data=[ecogHold.data ecog.data];
%                 ecogHold2.data=cat(3,ecogHold2.data, ecogSTFT.data);

                DataAfterCAR(:,:,finalPos+1:finalPos+size(LogNewData,3))=LogNewData;
                %DataAfterCAR2=cat(3,DataAfterCAR2, LogNewData);
                ANNewData_finalMAT(:,finalPosA+1:finalPosA+size(ANNewDataMAT,2))=ANNewDataMAT(:,1:size(ANNewDataMAT,2));
                %d=[ d ecog.data];
                if finalPos==0
                    specTimes(finalPos+1:finalPos+size(LogNewData,3))=ecogSTFT.time;
%                     specTimes2=ecogSTFT.time;
                else
                    specTimes(finalPos+1:finalPos+size(LogNewData,3))=ecogSTFT.time+specTimes(end);
%                     specTimes2=[specTimes2 specTimes2(end)+ecogSTFT.time];
                end
                finalPos=finalPos+size(LogNewData,3);
                finalPosA=finalPosA+size(ANNewDataMAT,2);
                ii=ii+1;
                %find events
                event=(ANNewData_finalMAT(desired_ANchan,:)>threshold);
                trigger=(diff(event)>0);
                detected_num_events = sum(trigger);
                handles=plotEventCounter(handles,'update',matlab_loopCounter,detected_num_events)
                
                %%
                if detected_num_events>0
                    %%% PARSE EVENTS HERE
                    [num_events,~,eventTimes] = ...
                        parseEvents(ANNewData_finalMAT, finalPos, desired_ANchan,...
                        desired_ANchan2, threshold, threshold2, analog.sampFreq, prev_num_events,...
                        [], [], number_of_analog,window_around_event/2);
                    
                    if num_events>prev_num_events
                        % num_events is #events after removing ones within 1 second, intialized to 0
                        prev_num_events = num_events;
                        newEventFlag=1;
                    else
                        newEventFlag=0;
                    end
                    
                    if  newEventFlag==1 && number_of_analog ==1
                        [eventRelatedAvg,num_samples,last_used_ind,allStacked,new_plot_flag,good_event_count] = average_event_window_angela(num_events, eventTimes(1,:), window_around_event,...
                            DataAfterCAR,finalPos(1),number_of_electrodes_total,...
                            num_avgEvent_freqs2plot, avgEvent_freqs2plot,...
                            eventRelatedAvg, num_samples, last_used_ind,allStacked,good_event_count, averages, stdevs,freq_band_singlestacked,specTimes);
                    elseif newEventFlag==1 && number_of_analog ==2
                        [eventRelatedAvg,last_used_ind,allStacked,new_plot_flag,current_num_events,lastSpec,lastSpec_event2,eventRelatedAvg2] =...
                            average_event_window_2_events(num_events, event_indices, window_around_event,window_after_event,...
                            DataAfterCAR,finalPos(1),number_electrodes ,num_freq_bands, freqs2plot,...
                            old_average, last_used_ind,allStacked,good_event_count, averages, stdevs,old_average_event2,freq_band_singlestacked)
                    end
                    
                    if  new_plot_flag==1 & newEventFlag==1
                        handles=plotEventCounter(handles,'newEvent',matlab_loopCounter,good_event_count)
                        %%
                        %Plot average Spectrogram(s) around event(s)
                        %always plot average spectrogram aligned
                        %around event 1
                        
                        if ave_spec_flag==1
                            %ZscoreEventRelatedAvg=(eventRelatedAvg-averages)./stdevs;
                            to_plot=allStacked(:,1:good_event_count,:,:);
                            to_plot=squeeze(nanmean(to_plot,2));
                            flattened = reshape_3Ddata(to_plot, window_around_event,size(to_plot,2), to_plot_grid);
                            to_plot=flattened;
                            handles=plotAveSpec(handles,'update',window_around_event,to_plot_grid,to_plot,number_of_electrodes_total,num_events,1)
                            if number_of_analog==2
                                ZscoreEventRelatedAvg=(eventRelatedAvg2-averages)./stdevs;
                                flattened = reshape_3Ddata(to_plot, window_around_event,size(ZscoreEventRelatedAvg,2), to_plot_grid);
                                handles=plotAveSpec(handles,flag,window_around_event,to_plot_grid,to_plot,number_of_electrodes_total,num_events,2)
                            end
                        end
                        
                        %plot single stacked
                        if single_stacked_flag ==1
                            to_plot=mean(allStacked(:,1:good_event_count,desired_freq_plot,:),3);
                            flattened = reshape_3Ddata(to_plot, window_around_event,good_event_count, to_plot_grid);
                            to_plot=flattened;
                            handles=plotStacked(handles,'update',window_around_event,to_plot_grid,to_plot,number_of_electrodes_total,good_event_count)
                        end
                        
                        %% plot mr
                        if ~isempty(get(handles.inputMRIFile,'String'))
                            %ZscoreEventRelatedAvg=(eventRelatedAvg-averages)./stdevs;
                            to_plot=allStacked(:,1:good_event_count,:,:);
                            to_plot=squeeze(nanmean(to_plot,2));
                            %flattened = reshape_3Ddata(to_plot, window_around_event,size(to_plot,2), to_plot_grid);
                            %to_plot=flattened;
                            handles=plotMRI(handles,'update',to_plot,[],sampPerSecond,desired_freq_plot);
                        end
                    end
                    
                end
                try
                    ANoldAcqPos=ANAcqPos;
                end
                matlab_loopCounter=matlab_loopCounter+1;
                try
                    oldAcqPos=AcqPos;
                end
                
                if toc< .1%integration_time%*sampling_rate
                    pause(.1)
                end
                tmp (ii)= toc;
            end
        else
            msgbox('Incorrect OpenEx Project')
        end
    else
        msgbox('OpenEx project Failed To Run')
    end
    %    end
else
    msgbox('OpenEx project not loaded reload OpenEx project and restart MATLAB script')
end
DA.CloseConnection
%
% disp('Saving figures...')
% if exist([subj_id '_avg.jpg'],'file')
%     subj_id = inputdlg('ID already exists, type new ID','',1,{subj_id});
%     if isempty(subj_id)
%         subj_id = {'EC_b'};
%     end
%     subj_id = subj_id{1};
% end
% saveas(fh,[subj_id '_avg.jpg'],'jpg')
% saveas(fh6,[subj_id '_stacked.jpg'],'jpg')
%
% play_movie = 0; %input('Play movie (1/0): ');
% while play_movie
%     save_movie = 0;
%     makeMovie(subj_id, ZscoreEventRelatedAvg,to_plot_grid, number_of_electrodes_total, save_movie);
%     play_movie = input('Play movie (1/0): ');
% end

disp('Done')


function AcqPos = AcquirePosition(sproj_name, BufferSize, DA)

for i = 1:length(sproj_name)
    AcqPos(i)=DA.GetTargetVal(sproj_name{i});%DA.GetTargetVal(sproj_name{i})-number_of_points2read+BufferSize(i);
    %AcqPos(i) = mod(AcqPos(i),BufferSize(i)); %%
    % + BufferSize takes care of circular buffer, prevents tdt from getting
    % confused over a negative position ( -number_of_points2read)
end

%%%%%%
function newbufferCounter = updateCounter(bufferCounter, AcqPos, oldAcqPos)
for i=1:length(AcqPos)
    newbufferCounter(i)=bufferCounter(i)+(AcqPos(i)<oldAcqPos(i));
end

%%%%%%
function [posInNewData] = findPosition(bufferCounter, BufferSize, AcqPos, number_of_channels)

for i=1:length(BufferSize)
    posInNewData(i) = (bufferCounter(i)*BufferSize(i)+AcqPos(i)-BufferSize(i))/number_of_channels;
end

%%%%%%
function  NewData=readTarget(dproj_name,AcqPos, number_of_points2read, DA)
for i = 1:length(AcqPos)
    NewData(i,:)=DA.ReadTargetVEX(dproj_name{i},AcqPos(1), number_of_points2read,'F32','F64');
end



function [average,num_samples,last_used_ind] = average_event_window(num_events, ind, window_around_event,...
    DataAfterCAR,AmountOfData,number_electrodes ,num_freq_bands, freqs2plot,...
    old_average, old_num_samples, last_used_ind)

%allocate memory, set up counter
average = zeros(number_electrodes,num_freq_bands, window_around_event+1);
sum_windows = zeros(number_electrodes,num_freq_bands, window_around_event+1);
num_samples = 0;

if isempty(old_average)
    old_average = zeros(number_electrodes,num_freq_bands, window_around_event+1);
end
if isempty(find(ind==last_used_ind))==1
    x=0;
else
    x=find(ind==last_used_ind);
end
%loop over each event starting from one after the last index, grab corresponding window
for i=(x+1):num_events
    beginning = ind(i) - window_around_event/2;
    last = ind(i)+ window_around_event/2;
    if beginning <1 || last>AmountOfData
        continue; %since can't add different sized window, just ignore
    end
    timepts = beginning:last;
    sum_windows = sum_windows + DataAfterCAR(:,freqs2plot, timepts);
    num_samples = num_samples + 1;
    last_used_ind = ind(i);
end;
if num_samples>0
    num_samples = num_samples+old_num_samples;
    average = (old_average*old_num_samples + sum_windows)/(num_samples);
else
    average = old_average;
    num_samples = old_num_samples;
end
%%%%




function single_window = single_event_window(ind, window_around_event,...
    DataAfterCAR,AmountOfData,number_electrodes,num_freq_bands, freqs2plot)

single_window=zeros(number_electrodes,num_freq_bands, window_around_event+1);
if length(ind)>1
    for i=0:1
        beginning = ind(end-i)-window_around_event/2;
        last = ind(end-i)+window_around_event/2;
        if beginning <1 || last > AmountOfData
            continue %since can't add different sized window, just ignore
        else
            timepts = beginning:last;
            single_window=DataAfterCAR(:,freqs2plot,timepts);
            return
        end
    end
end





function to_plot = reshape_3Ddata(to_reshape, window_around_event,...
    num_freq_bands, to_plot_grid)
%to_reshape = chan x freq x timepts
to_plot = zeros(num_freq_bands*to_plot_grid, to_plot_grid*(window_around_event+1));

beginning = 0;
for i = 0:(to_plot_grid-1)  %ASSUMES to_plot_grid x to_plot_grid CHANNELS
    last =  beginning+to_plot_grid*(window_around_event+1);
    data = to_reshape(i*to_plot_grid+1:i*to_plot_grid+to_plot_grid,:,:);
    to_plot(i*num_freq_bands+1:i*num_freq_bands+num_freq_bands,:) = ...
        flipud(reshape(shiftdim(data,1),num_freq_bands, to_plot_grid*(window_around_event+1)));
    beginning = last;
    %Grab every eight channels and form to_plot_grid x to_plot_grid square
end



function [baseline,averages, stdevs, medians,logBaselineDataMAT]=calculate_baseline(time2collectBaseline, ...
    sproj_name,dproj_name, points_needed4baseline, BufferSize, DA, num_channels, num_freq_bands,...
    sampling_rate, number_of_sec2read)

h = waitbar(0,'Calculating Baseline...');
for i = 1:time2collectBaseline
    pause(1);
    waitbar(i/(time2collectBaseline+1),h)
end
pause(1);
AcqPos=AcquirePosition(sproj_name, BufferSize, DA);
AcqPos = AcqPos - points_needed4baseline+BufferSize(1);
AcqPos = mod(AcqPos,BufferSize(1));

BaselineData =readTarget(dproj_name,AcqPos, points_needed4baseline, DA);
% converts to matrix
% reshape to chan x freq x time
BaselineDataMAT=[];
num_freq_bands=1;
for i = 1:size(BaselineData,1)
    BaselineDataMAT= [BaselineDataMAT; shiftdim(reshape(reshape(BaselineData(i,:),...
        num_freq_bands*num_channels, sampling_rate*time2collectBaseline)',...
        sampling_rate*time2collectBaseline,num_channels,num_freq_bands),1)];
end

medians = median(BaselineDataMAT,3);
%logBaselineDataMAT = log(BaselineDataMAT+repmat(medians,[1 1 size(BaselineDataMAT,3)])+eps);

stdevs=[];
medians=[];
averages=[];
logBaselineDataMAT=[];
%medians=repmat(medians,[1 1 (sampling_rate*number_of_sec2read)]);
%averages = repmat(mean(logBaselineDataMAT,3),[1 1 (sampling_rate*number_of_sec2read+1)]);
%stdevs = repmat(std(logBaselineDataMAT,0,3),[1 1 (sampling_rate*number_of_sec2read+1)]);
%averages = repmat(mean(logBaselineDataMAT,3),[1 1 (sampling_rate*1+1)]);
%stdevs = repmat(std(logBaselineDataMAT,0,3),[1 1 (sampling_rate*1+1)]);
waitbar(1,h,'Baseline Statistics Gathered!')
pause(1)
close(h)
baseline.data=BaselineDataMAT;
baseline.sampFreq=sampling_rate;

function to_plot = SingleStackedTrials(DataMat, ind, num_events,...
    averages, stdevs, freq_band, window_around_event)
% function to_plot = SingleStackedTrials(DataMat, ind, ...
%   averages, stdevs, freq_band, window_around_event)
% Creates z-scored matrix (num_events x time) for ONE channel. Each row
% represents zscore data for one event
%
% DataMat - Chan x freq x time (envelope of raw data)
% AnalogData - event data
% AnalogDataChan
% freq_band
% windowaroundevent
% averages, stdevs = baseline stats

% format data, take log of it, extract Chan, freq_band
Data = DataMat(:, freq_band,:);

% format statistics, extract Chan, freq_band
averages = squeeze(averages(:, freq_band, :));
stdevs = squeeze(stdevs(:, freq_band,:));


%allocate memory
AmountOfData = length(Data);
to_plot = zeros(size(Data,1),num_events,window_around_event+1);
empty_trials = [];

%loop through every event and extract data around event
for i=1:num_events
    beginning = ind(i) - window_around_event/2;
    last = ind(i)+ window_around_event/2;
    if beginning <1 || last>AmountOfData
        empty_trials = [empty_trials i];
        continue; %since can't use different sized window, just ignore
    end
    timepts = beginning:last;
    zscored_logGdat = (Data(:,timepts)-averages)./stdevs;
    to_plot(:,i,:) = zscored_logGdat;
    
end

to_plot(:,empty_trials, :) = [];


function lastSpec = artifactrejection2(lastSpec, window_around_event,...
    min_zscore, max_zscore, percentage_min, percentage_max, to_plot_grid)
%reshaped data = freqs2plot*8 x window_around_event+1*8
%
%  if more than percentage_min% of the data is below min_zscore, reject
%  if more than percentage_max%of the data is above max_zscore, reject

[a,b]=find(abs(mean(lastSpec,3))>max_zscore);
for i=1:length(a)
    lastSpec(a(i),b(i),:)=NaN;
end

function to_plot = artifactrejection(reshaped_data, window_around_event,...
    min_zscore, max_zscore, percentage_min, percentage_max, to_plot_grid)
%reshaped data = freqs2plot*8 x window_around_event+1*8
%
%  if more than percentage_min% of the data is below min_zscore, reject
%  if more than percentage_max%of the data is above max_zscore, reject

rows = size(reshaped_data,1);
cols = 1:window_around_event+1:to_plot_grid*window_around_event;
to_plot = reshaped_data;
for i = 1:rows
    for j = 1:to_plot_grid
        k = cols(j):cols(j)+(window_around_event+1)-1;
        if (sum(reshaped_data(i,k)<min_zscore)>=percentage_min*(window_around_event+1) ||...
                sum(reshaped_data(i,k)>max_zscore)>=percentage_max*(window_around_event+1))
            to_plot(i,k)=zeros(1,window_around_event+1);
        end
    end
end


% --- Executes on selection change in get_total_electrodes.
function get_total_electrodes_Callback(hObject, eventdata, handles)
% hObject    handle to get_total_electrodes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns get_total_electrodes contents as cell array
%        contents{get(hObject,'Value')} returns selected item from get_total_electrodes
handles.number_of_electrodes_total = cellstr(get(hObject,'String'))
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function get_total_electrodes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to get_total_electrodes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%
function time2collectBaseline_Callback(hObject, eventdata, handles)
% hObject    handle to time2collectBaseline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of time2collectBaseline as text
%        str2double(get(hObject,'String')) returns contents of time2collectBaseline as a double


% --- Executes during object creation, after setting all properties.
function time2collectBaseline_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time2collectBaseline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function desired_freq_plot_Callback(hObject, eventdata, handles)
% hObject    handle to desired_freq_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of desired_freq_plot as text
%        str2double(get(hObject,'String')) returns contents of desired_freq_plot as a double


% --- Executes during object creation, after setting all properties.
function desired_freq_plot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to desired_freq_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in calculated_baseline_flag.
function calculated_baseline_flag_Callback(hObject, eventdata, handles)
% hObject    handle to calculated_baseline_flag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of calculated_baseline_flag



function threshold_Callback(hObject, eventdata, handles)
% hObject    handle to threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of threshold as text
%        str2double(get(hObject,'String')) returns contents of threshold as a double


% --- Executes during object creation, after setting all properties.
function threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function freq_band_singlestacked_Callback(hObject, eventdata, handles)
% hObject    handle to freq_band_singlestacked (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of freq_band_singlestacked as text
%        str2double(get(hObject,'String')) returns contents of freq_band_singlestacked as a double


% --- Executes during object creation, after setting all properties.
function freq_band_singlestacked_CreateFcn(hObject, eventdata, handles)
% hObject    handle to freq_band_singlestacked (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function subj_id_Callback(hObject, eventdata, handles)
% hObject    handle to subj_id (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of subj_id as text
%        str2double(get(hObject,'String')) returns contents of subj_id as a double


% --- Executes during object creation, after setting all properties.
function subj_id_CreateFcn(hObject, eventdata, handles)
% hObject    handle to subj_id (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function desired_ANchan_Callback(hObject, eventdata, handles)
% hObject    handle to desired_ANchan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of desired_ANchan as text
%        str2double(get(hObject,'String')) returns contents of desired_ANchan as a double


% --- Executes during object creation, after setting all properties.
function desired_ANchan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to desired_ANchan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in numchannels.
function numchannels_Callback(hObject, eventdata, handles)
% hObject    handle to numchannels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns numchannels contents as cell array
%        contents{get(hObject,'Value')} returns selected item from numchannels


% --- Executes during object creation, after setting all properties.
function numchannels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numchannels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in number_of_electrodes_total.
function number_of_electrodes_total_Callback(hObject, eventdata, handles)
% hObject    handle to number_of_electrodes_total (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns number_of_electrodes_total contents as cell array
%        contents{get(hObject,'Value')} returns selected item from number_of_electrodes_total


% --- Executes during object creation, after setting all properties.
function number_of_electrodes_total_CreateFcn(hObject, eventdata, handles)
% hObject    handle to number_of_electrodes_total (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over number_of_electrodes_total.
function number_of_electrodes_total_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to number_of_electrodes_total (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on update and none of its controls.
function update_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to update (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on mouse press over axes background.
function eventCounter_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to eventCounter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in stop_flag.
function stop_flag_Callback(hObject, eventdata, handles)
% hObject    handle to stop_flag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.DA.SetSysMode(1)% 1 = Standby, 2 = Preview


% --- Executes on button press in mri_flag.
function mri_flag_Callback(hObject, eventdata, handles)
% hObject    handle to mri_flag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of mri_flag


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    load(['C:\Users\Angela_2\Dropbox\Connie_Angela (1)\EC9_B58\TDT_DataMAT' '\EC9_B58_start_ind_AND_threshold.mat'])
    load(['C:\Users\Angela_2\Dropbox\Connie_Angela (1)\EC9_B58\TDT_DataMAT' '\DataMAT.mat'])
    load(['C:\Users\Angela_2\Dropbox\Connie_Angela (1)\EC9_B58\TDT_DataMAT' '\baseline_stats_EC9_b58.mat'])
catch
    load(['C:\Users\ego\Dropbox\Connie_Angela (1)\EC9_B58\TDT_DataMAT' '\EC9_B58_start_ind_AND_threshold.mat'])
    load(['C:\Users\ego\Dropbox\Connie_Angela (1)\EC9_B58\TDT_DataMAT' '\DataMAT.mat'])
    load(['C:\Users\ego\Dropbox\Connie_Angela (1)\EC9_B58\TDT_DataMAT' '\baseline_stats_EC9_b58.mat'])
end
clear threshold

vars=whos;
for i=1:length(vars)
    handles.(vars(i).name)=eval(vars(i).name);
end
handles.simData=1;
fprintf('Done\n')
guidata(hObject, handles);


function  [LogNewData,ANNewData,index,logGdat,analog_dat,ANAcqPos,AcqPos,bufferCounter]=getSimData(handles,time2collectBaseline,sampling_rate,number_of_sec2read,ii,matlab_loopCounter,bufferCounter);
points_needed4baseline = sampling_rate*time2collectBaseline;
number_of_points2read = sampling_rate*number_of_sec2read;

Envelope=handles.Envelope;
AnalogData=handles.AnalogData;
medians=handles.medians;
start_ind=handles.start_ind;
f=fields(handles);
averages=handles.averages;
stdevs=handles.stdevs;

gdat = Envelope(:,:,start_ind:end);
analog_dat = AnalogData(:,start_ind:end);
index = 1:number_of_points2read:length(gdat);

clear Envelope AnalogData RawData
m = medians(:,:,1);
m = repmat(m, [1 1 length(gdat)]);
logGdat = log(abs(gdat)+m+eps);
LogNewData = logGdat(:,:,index(ii):index(ii+1)-1);
ANNewData= analog_dat(:,1:index(ii+1)-1);
DataAfterCAR = logGdat(:,:,1:index(ii+1)-1);
ANfinalPos=index(ii+1)-1;
ANAcqPos=0;
AcqPos=0;


function  [ecog,analog]=loadHTKData(handles,number_of_sec2read,ii)
filenames=getFileNames('blocks',1:256);
ecog=loadHTKFile_Name([get(handles.htkFilePath,'String') filesep 'Downsampled400'],filenames);
%ecog=downsampleEcog(ecog,400,ecog.sampFreq);
analog=loadHTKFile_Name([get(handles.htkFilePath,'String') filesep 'Analog'],{'ANIN1','ANIN2','ANIN3','ANIN4'});

for i=1:size(analog.data,1)
    analog.data(i,:)=abs(hilbert(analog.data(i,:)));
end
analog=downsampleEcog(analog,ecog.sampFreq,analog.sampFreq);


function [ecog,analog,currentTime]=grabDataBuffer(ecogAll,analogAll,lastTime)
if lastTime==0
    elapsedTime=2;
else
    elapsedTime=toc;
end
tic;
samps=[round(lastTime*ecogAll.sampFreq)+1:round((lastTime+elapsedTime)*ecogAll.sampFreq)];
ecog.data=ecogAll.data(:,samps);
ecog.sampFreq=ecogAll.sampFreq;
analog.data=analogAll.data(:,samps);
analog.sampFreq=analogAll.sampFreq;
currentTime=lastTime+elapsedTime;


function ecog=loadHTKBaseline(handles,time2collectBaseline)
filenames=getFileNames('blocks',1:256);
ecog=loadHTKFile_Name([get(handles.inputBaselineFile,'String') filesep 'Downsampled400'],filenames,[0 time2collectBaseline]*1000);
%ecog=downsampleEcog(ecog,400,ecog.sampFreq);

function  [ecog,analog,LogNewData,ANNewData,index,logGdat,analog_dat,ANAcqPos,AcqPos,bufferCounter, ...
    constant, ANconstant, posInNewData_AfterCAR, ANposInNewData,number_of_points2read,ANnumber_of_points2read]=getTDTData(handles,sampling_rate,ii,matlab_loopCounter,bufferCounter,ANbufferCounter,...
    number_of_points2read, ANnumber_of_points2read,sproj_name,dproj_name,sANproj_name,dANproj_name,...
    BufferSize,ANBufferSize, oldAcqPos, ANoldAcqPos,number_of_channels,ANnumber_of_channels, constant,ANconstant,...
    posInNewData_AfterCAR, ANposInNewData)
LogNewData = [];index=[];logGdat=[];analog_dat=[];
DA=handles.DA;
if matlab_loopCounter ~= 0
    AcqPos = oldAcqPos+number_of_points2read;
    AcqPos = mod(AcqPos,BufferSize(1));
    tmp=AcquirePosition(sproj_name, BufferSize, DA);
    number_of_points2read = tmp(1)-AcqPos(1);
    if number_of_points2read<0
        number_of_points2read = BufferSize(1) - AcqPos(1)+tmp(1);
    end
else
    AcqPos=AcquirePosition(sproj_name, BufferSize, DA) - number_of_points2read +BufferSize(1);
    AcqPos = mod(AcqPos,BufferSize(1));
end
num_freq_bands = 1;
bufferCounter=updateCounter(bufferCounter, AcqPos, oldAcqPos);
posInNewData_AfterCAR=findPosition(bufferCounter,BufferSize, AcqPos, ...
    number_of_channels*num_freq_bands);
if matlab_loopCounter~=0
    posInNewData_AfterCAR=posInNewData_AfterCAR-constant;
else
    % first position may be late in the buffer,
    % this sets intial matlab buffer index to one
    constant = posInNewData_AfterCAR-1;
    posInNewData_AfterCAR = ones(size(posInNewData_AfterCAR));
end

% Repeat for Analog channels
    if matlab_loopCounter ~= 0
        ANAcqPos = ANoldAcqPos+ANnumber_of_points2read;
        ANAcqPos = mod(ANAcqPos,ANBufferSize);
        tmp=AcquirePosition(sANproj_name, ANBufferSize, DA);
        ANnumber_of_points2read = tmp(1)-ANAcqPos(1);
        if ANnumber_of_points2read<0
            ANnumber_of_points2read = ANBufferSize - ANAcqPos(1)+tmp(1);
        end
    else
        ANAcqPos=AcquirePosition(sANproj_name, ANBufferSize, DA) - ANnumber_of_points2read+ ANBufferSize;
        ANAcqPos = mod(ANAcqPos,ANBufferSize);
    end
    
    ANbufferCounter = updateCounter(ANbufferCounter, ANAcqPos, ANoldAcqPos);
    ANposInNewData=findPosition(ANbufferCounter, ANBufferSize, ANAcqPos, ANnumber_of_channels);
    if matlab_loopCounter~=0
        ANposInNewData=ANposInNewData-ANconstant;
    else
        ANconstant = ANposInNewData-1;
        ANposInNewData = ones(size(ANposInNewData));
    end
    ANNewData = DA.ReadTargetVEX(dANproj_name{1},ANAcqPos, ANnumber_of_points2read,'F32','F64');
    ANNewDataMAT=reshape(ANNewData,ANnumber_of_channels,ANnumber_of_points2read/ANnumber_of_channels);

% Read from TDT buffers and reshape to chan x freq bands x time
NewData=readTarget(dproj_name, AcqPos, number_of_points2read, DA);
NewDataMAT=zeros(number_of_channels*length(AcqPos),num_freq_bands,number_of_points2read/(num_freq_bands*number_of_channels));
for i = 1:size(NewData,1)
    ind = (i-1)*number_of_channels+1;
    NewDataMAT( ind:(ind+number_of_channels-1),:, :) =shiftdim(reshape(reshape(NewData(i,:),...
        num_freq_bands*number_of_channels, number_of_points2read/(num_freq_bands*number_of_channels))',...
        number_of_points2read/(num_freq_bands*number_of_channels),number_of_channels,num_freq_bands),1);
end

ecog.data=squeeze(NewDataMAT);
ecog.sampFreq=sampling_rate;

analog.data=ANNewDataMAT;
analog.sampFreq=sampling_rate;

%Take log(power+medians)
%tmp = repmat(medians,[1 1 (number_of_points2read/(num_freq_bands*number_of_channels))]);
%LogNewData =log(abs(NewDataMAT)+tmp+eps);

function handles=plotMRI(handles,flag,ZscoreEventRelatedAvg,filePath,sampPerSecond,hgfreq)

switch flag
    case 'initiate'
        % MRI_FLAG PLOT
        set(0,'currentfigure',handles.figure1);
        set(handles.figure1,'CurrentAxes',handles.mriPlot1); cla
        set(handles.figure1,'CurrentAxes',handles.mriPlot2); cla
        set(handles.figure1,'CurrentAxes',handles.mriPlot3); cla
        set(handles.figure1,'CurrentAxes',handles.mriPlot4); cla
        fhmr{1}=handles.mriPlot1; set(fhmr{1},'XTickLabel',[]); set(fhmr{1},'YTickLabel',[])
        fhmr{2}=handles.mriPlot2; set(fhmr{2},'XTickLabel',[]); set(fhmr{2},'YTickLabel',[])
        fhmr{3}=handles.mriPlot3; set(fhmr{3},'XTickLabel',[]); set(fhmr{3},'YTickLabel',[])
        fhmr{4}=handles.mriPlot4; set(fhmr{4},'XTickLabel',[]); set(fhmr{4},'YTickLabel',[])
        colormapjet = colormap('jet');
        
        
        %tp_toplot = {-500:2:-250; -250:2:0; 0:2:250; 250:2:500};
        
        %load regdata
        %[fname, pth] = uigetfile('*.jpg','Load MRI jpg');
        img = imread(filePath);
        [pth, fname, tmp] = fileparts(filePath);
        fname = [pth filesep 'regdata.mat'];
        if exist(fname,'file')
            regdata = load(fname);
        else
            [fname, pth ] = uigetfile('*regdata.mat', 'Load electrode positions');
            regdata = load([pth fname]);
        end
        
        set(0,'currentfigure',handles.figure1);
        for j = 1:4%length(tp_toplot)
            %tp_toplot{j} =  tp_toplot{j}/1000*500 + 251;
            
            set(handles.figure1,'CurrentAxes',fhmr{j});
            imagesc(img); hold on
            
            mrihandle{j} = scatter(regdata.xy(1,:),regdata.xy(2,:),15,'filled');
            set(fhmr{j},'XTickLabel',[]); set(fhmr{j},'YTickLabel',[])
        end
        %handles.MRIparams.tp_toplot=tp_toplot;
        handles.MRIparams.fhmr=fhmr;
        %handles.MRIparams.hgfreq=hgfreq;
        handles.MRIparams.colormapjet=colormapjet;
        handles.MRIparams.regdata=regdata;
    case 'update'
        set(0,'currentfigure',handles.figure1);
        f=fieldnames(handles.MRIparams);
        for i=1:length(f)
            eval(sprintf('%s=deal(handles.MRIparams.(f{i}));',f{i}));
        end
        %hgfreq = 8:13;
        sampInterval=round(size(ZscoreEventRelatedAvg,3)/2)-round(sampPerSecond/2):round(sampPerSecond/4):...
            round(size(ZscoreEventRelatedAvg,3)/2)+round(sampPerSecond/2);
        tp_toplot={sampInterval(1):sampInterval(2);sampInterval(2):sampInterval(3);...
            sampInterval(3):sampInterval(4);sampInterval(4):sampInterval(5)};
        for j = 1:length(tp_toplot)
            try
                set(handles.figure1,'CurrentAxes',fhmr{j});
                to_plot = ZscoreEventRelatedAvg(:,hgfreq,tp_toplot{j});
                mr1data = mean(mean(to_plot,3),2);
                cmax = 3; cmin = -3;
                %cmax=max(mr1data); cmin=min(mr1data);
                colormap_index = fix((mr1data-cmin)/(cmax-cmin)*length(colormapjet))+1;
                colormap_index(colormap_index>length(colormapjet)) = length(colormapjet);
                colormap_index(colormap_index<=0) = 1;
                try
                    set(mrihandle{j},'CData',colormapjet(colormap_index,:))
                catch
                    mrihandle{j} = scatter(regdata.xy(1,:),regdata.xy(2,:),15,colormapjet(colormap_index,:),'filled');
                    handles.MRIparams.mrihandle=mrihandle;
                end
                
                set(fhmr{j},'XTickLabel',[]); set(fhmr{j},'YTickLabel',[])
            end
        end
end

function handles=plotStacked(handles,flag,window_around_event,to_plot_grid,to_plot,number_of_electrodes_total,num_events)
pcolorPlot=1;
switch flag
    case 'initiate'
        %SINGLE STACKED PLOT
        f=fieldnames(handles.numberingParams);
        for i=1:length(f)
            eval(sprintf('%s=deal(handles.numberingParams.(f{i}));',f{i}));
        end
        fh6 = figure(6); clf %single stacked plot
        z=axes('Position',[.1 .1 .85 .85],'visible','off');
        set(fh6,'DefaultAxesxtick',[0.5: (window_around_event+1):to_plot_grid*(window_around_event+1)])
        set(fh6,'DefaultAxeslinewidth',1)
        set(fh6,'DefaultAxesyticklabel','')
        set(fh6,'DefaultAxesxticklabel','')
        set(fh6,'DefaultAxesgridlinestyle','-')
        set(fh6, 'Name','Single Stacked','NumberTitle','off')
        set(fh6,'CurrentAxes',z);
        handles.StackedParams.fh6=fh6;
        set(0,'currentfigure',fh6);
    case 'update'
        set(0,'currentfigure',handles.StackedParams.fh6);
        f=fieldnames(handles.numberingParams);
        for i=1:length(f)
            eval(sprintf('%s=deal(handles.numberingParams.(f{i}));',f{i}));
        end
        if pcolorPlot==0
            imagesc(to_plot,[-7 7])
            if num_events >1
                fig_num_y_stacked = [1.1:num_events:to_plot_grid*num_events];
                fig_num_y_stacked = reshape(repmat(fig_num_y_stacked, [to_plot_grid 1]), number_of_electrodes_total, 1);
                stackedfig_dash_y = repmat([0 num_events*to_plot_grid],to_plot_grid,1)';
                line(fig_dash_x,stackedfig_dash_y,'LineStyle','--','linewidth',1','color',[0 0 0]);
                set(gca,'YTick',[0.5: num_events:to_plot_grid*num_events])
                text(fig_num_x, fig_num_y_stacked, fig_nums)
            end
        else
            pcolor(flipud(to_plot));
            shading interp;
            set(gca,'CLim',[-7 7]);
            if num_events >1
                fig_num_y_stacked = [1.1:num_events:to_plot_grid*num_events];
                fig_num_y_stacked = reshape(repmat(fig_num_y_stacked, [to_plot_grid 1]), number_of_electrodes_total, 1);
                stackedfig_dash_y = repmat([0 num_events*to_plot_grid],to_plot_grid,1)';
                line(fig_dash_x,stackedfig_dash_y,'LineStyle','--','linewidth',1','color',[0 0 0]);
                set(gca,'YTick',[0.5: num_events:to_plot_grid*num_events])
                text(fig_num_x, flipud(fig_num_y_stacked),repmat(1,size(fig_num_y_stacked)), fig_nums);
            end
        end
        set(gca,'layer','top');
        grid on
        drawnow
end


function handles=plotContinual(handles,flag,window_around_event,to_plot_grid,to_plot,number_of_electrodes_total,num_events)
switch flag
    case 'initiate'
        set(0,'currentfigure',handles.figure1);
        set(handles.figure1,'CurrentAxes',handles.continualAverage); cla
        fh2=handles.continualAverage;
        title('Continual Plot')
        set(fh2,'XTickLabel',[]);
        set(fh2,'YTickLabel',[]);
        handles.ContinualPlot.fh2=fh2;
    case 'update'
        fh2=handles.ContinualPlot.fh2;
        f=fieldnames(handles.numberingParams);
        for i=1:length(f)
            eval(sprintf('%s=deal(handles.numberingParams.(f{i}));',f{i}));
        end
        set(0,'currentfigure',handles.figure1);
        set(handles.figure1,'CurrentAxes',fh2)
        try
            set(handles.ContinualPlot.imageH,'CData',to_plot);
        catch
            handles.ContinualPlot.imageH=imagesc(to_plot, [-7 7]);
        end
        %imagesc(to_plot, [-7 7]);
        drawnow;
end


function handles=plotAveSpec(handles,flag,window_around_event,to_plot_grid,to_plot,number_of_electrodes_total,num_events,number_of_analog,num_avgEvent_freqs2plot)
pcolorPlot=1;

switch flag
    case 'initiate'
        %SINGLE STACKED PLOT
        f=fieldnames(handles.numberingParams);
        for i=1:length(f)
            eval(sprintf('%s=deal(handles.numberingParams.(f{i}));',f{i}));
        end
        for i=1:number_of_analog
            fh=figure(i+1); clf% figure for average event trigger analysis
            z=axes('Position',[.1 .1 .85 .85],'visible','off');
            set(fh,'DefaultAxesxtick',[0.5: (window_around_event+1):to_plot_grid*(window_around_event+1)])
            set(fh,'DefaultAxesytick',[0.5: num_avgEvent_freqs2plot:to_plot_grid*num_avgEvent_freqs2plot])
            set(fh,'DefaultAxesyticklabel','')
            set(fh,'DefaultAxesxticklabel','')
            set(fh,'DefaultAxeslinewidth',1)
            set(fh,'DefaultAxesgridlinestyle','-')
            set(fh, 'Name','Average Event Time Frequency Plot','NumberTitle','off')
            set(fh,'CurrentAxes',z);
            handles.aveSpec(i).fh=fh;
            
            if pcolorPlot==0
                handles.aveSpec(i).imhandle=imagesc(to_plot,[-3 3]);
                text(fig_num_x, fig_num_y, fig_nums)
                line(fig_dash_x,avgfig_dash_y,'LineStyle','--','linewidth',1','color',[0 0 0]);
            else
                handles.aveSpec(i).imhandle=pcolor(flipud(to_plot));
                shading interp
                set(gca,'CLim',[-3 3]);
                text(fig_num_x, flipud(fig_num_y), repmat(1,size(fig_num_x)),fig_nums)
                line(fig_dash_x,avgfig_dash_y,repmat(1,size(avgfig_dash_y)),'LineStyle','--','linewidth',1','color',[0 0 0]);
            end
            set(gca,'layer','top')
            grid on;
        end
    case 'update'
        %fh=handles.aveSpec(number_of_analog).fh;
        %set(0,'currentfigure',fh)
        if pcolorPlot==0
            set(handles.aveSpec(number_of_analog).imhandle,'Cdata',(to_plot));
        else
            set(handles.aveSpec(number_of_analog).imhandle,'Cdata',flipud(to_plot));
        end
        
        drawnow;
end


function handles=plotEventCounter(handles,flag,matlab_loopCounter,detected_num_events)
switch flag
    case 'initiate'
        %SINGLE STACKED PLOT
        set(0,'currentfigure',handles.figure1);
        set(handles.figure1,'CurrentAxes',handles.eventCounter); cla
        fh5=handles.eventCounter;
        set(fh5,'XTick',[0:100])
        set(fh5,'YTick',[0:100])
        set(fh5,'XTickLabel',[])
        handles.eventCounterParams.fh5=fh5;
    case 'update'
        set(0,'currentfigure',handles.figure1);
        f=fieldnames(handles.numberingParams);
        for i=1:length(f)
            eval(sprintf('%s=deal(handles.numberingParams.(f{i}));',f{i}));
        end
        fh5=handles.eventCounterParams.fh5;
        set(handles.figure1,'CurrentAxes',fh5); hold on
        plot(matlab_loopCounter,detected_num_events,'*');
        grid on;
    case 'newEvent'
        fh5=handles.eventCounterParams.fh5;
        f=fieldnames(handles.numberingParams);
        for i=1:length(f)
            eval(sprintf('%s=deal(handles.numberingParams.(f{i}));',f{i}));
        end
        set(0,'currentfigure',handles.figure1);
        set(handles.figure1,'CurrentAxes',fh5)
        hold on
        plot(matlab_loopCounter,detected_num_events,'r*');
end

function handles=getElectrodeNums(handles,num_avgEvent_freqs2plot,num_singleEvent_freqs2plot,number_of_electrodes_total,window_around_event,to_plot_grid)
% coordinates for the electrode numbers
fig_num_y = [1.1:num_avgEvent_freqs2plot:to_plot_grid*num_avgEvent_freqs2plot];
fig_num_y = reshape(repmat(fig_num_y,[to_plot_grid 1]),number_of_electrodes_total,1);
fig_num_x = [.5:(window_around_event+1):to_plot_grid*(window_around_event+1)];
fig_num_x = reshape(repmat(fig_num_x,[1 to_plot_grid]),number_of_electrodes_total,1);
fig_nums = {};
for i = 1:number_of_electrodes_total
    fig_nums = [fig_nums {num2str(i)}];
end

% coordinates for the dash lines
avgfig_dash_y = repmat([0 num_avgEvent_freqs2plot*to_plot_grid],to_plot_grid,1)';
sigfig_dash_y = repmat([0 num_singleEvent_freqs2plot*to_plot_grid],to_plot_grid,1)';
fig_dash_x = repmat([round((window_around_event+1)/2):window_around_event+1:to_plot_grid*(window_around_event+1)],2,1);


f=whos;
f={f.name};
f=f(find(~strcmp(f,'handles')));
for i=1:length(f)
    eval(['handles.numberingParams.(f{i})=' f{i}]);
end


function [average,num_samples,last_used_ind,allStacked,new_plot_flag,currentEvent,lastSpec] = average_event_window_angela(num_events, eventTimes, window_around_event,...
    DataAfterCAR,AmountOfData,number_electrodes ,num_freq_bands, freqs2plot,...
    old_average, old_num_samples, last_used_ind,allStacked,good_event_count, averages, stdevs,freq_band_singlestacked,specTimes)
lastSpec=[];
%allocate memory, set up counter
current_num_events=num_events;
% average = zeros(number_electrodes,num_freq_bands, window_around_event+1);
%sum_windows = zeros(number_electrodes,num_freq_bands, window_around_event+1);
num_samples = 0;
ind=findNearest(eventTimes,specTimes);
if isempty(old_average)
    old_average = NaN(number_electrodes,num_freq_bands, window_around_event+1);
end
if isempty(find(ind==last_used_ind))
    x=0;
else
    x=find(ind==last_used_ind);
end
new_plot_flag=0;
%loop over each event, grab corresponding window
currentEvent=x+1;
for  i=currentEvent:current_num_events
    beginning = ind(i) - floor(window_around_event/2);
    last = ind(i)+ ceil(window_around_event/2);
    if beginning <1 || last>AmountOfData
        continue; %since can't add different sized window, just ignore
    else
        timepts = beginning:last;
        lastSpec=(DataAfterCAR(:,freqs2plot, timepts)-averages)./stdevs;
        %sum_windows = sum_windows + DataAfterCAR(:,freqs2plot, timepts);
        num_samples = num_samples + 1;
        %last_used_ind = ind(currentEvent);
        lastSpec = artifactrejection2(lastSpec, window_around_event,-5.5, 5.5, .5, .8);
        allStacked(:,currentEvent,:,:)=lastSpec(:,:,:);
        new_plot_flag=1;
        currentEvent=currentEvent+1;
    end
end;
currentEvent=currentEvent-1;
% if num_samples>0
%     num_samples = num_samples+old_num_samples;
%     average = (old_average*old_num_samples + sum_windows)/(num_samples);
% else %if no new events
%     average = old_average;
%     num_samples = old_num_samples;
%
% end
average=[];


function [average,last_used_ind,allStacked,new_plot_flag, plotted_num_events,lastSpec,lastSpec_event2,average_event2,event_indices] =...
    average_event_window_2_events(num_events, event_indices, window_before_event,window_after_event,...
    DataAfterCAR,AmountOfData,number_electrodes ,num_freq_bands, freqs2plot,...
    old_average, last_used_ind,allStacked,averages, stdevs,old_average_event2,freq_band_singlestacked,plotted_num_events)
%allStacked is the single stacked variable where each trial is sorted by
%latency of response. Aligned at event1, and includes event2
%old_average_event2 is updated average spectrogram aligned at event2
lastSpec=[];
lastSpec_event2=[];
window_around_event=window_before_event+window_after_event;

ind=event_indices(1,:);
latency=event_indices(2,:)-ind;
%allocate memory, set up counter
%current_num_events=good_event_count;
sum_windows = zeros(number_electrodes,num_freq_bands, window_before_event*2+1);
sum_windows_event2 = zeros(number_electrodes,num_freq_bands, window_before_event*2+1);

if isempty(old_average)
    old_average = zeros(number_electrodes,num_freq_bands, window_before_event*2+1);
    old_average_event2=zeros(number_electrodes,num_freq_bands, window_before_event*2+1);
end


%Finds event segments within bounds of available data
startidx=plotted_num_events+1;
beginning = ind(startidx) - window_before_event;
endidx=num_events;
last = event_indices(2,endidx)+ window_before_event;
last2 = ind(endidx)+ window_after_event;
while((last>AmountOfData) | (last2>AmountOfData)) & endidx>=startidx% what if ind(num_events-1)+window/2 is still out of bounds?
    last = event_indices(2,endidx)+ window_before_event;
    last2 = ind(endidx)+ window_after_event;
    endidx=num_events-1;
end

%If there are new events within the bounds of current data
%do i even need this? will endidx always be >= startidx?
if endidx>=startidx
    %loop over each event, grab corresponding window
    for i=[startidx:endidx]
        beginning = ind(i) - window_before_event;
        last = ind(i)+ window_after_event; % NOTE: need different last for stacked, see ln 761.
        %I think we should grab 2.5seconds after the first event
        %new_count=new_count+1;
        timepts = beginning:last;
        lastSpec=DataAfterCAR(:,freqs2plot, timepts);
        last_used_ind = ind(i);
        %prev_num_events=prev_num_events+1;
        lastSpec=(lastSpec-averages)./stdevs;
        sum_windows = sum_windows + lastSpec(:,:,1:window_before_event*2+1);
        
        allStacked(:,i,:)=lastSpec(:,freq_band_singlestacked,:);
        
        %Get data aligned at 2nd event
        beginning = event_indices(2,i) - window_before_event;
        last = event_indices(2,i)+ window_before_event;
        timepts = beginning:last;
        
        lastSpec_event2=(DataAfterCAR(:,freqs2plot, timepts)-averages(:,:,1:2*window_before_event+1))./stdevs(:,:,1:2*window_before_event+1);
        sum_windows_event2 = sum_windows_event2 + lastSpec_event2;
        
    end;
    %sort stacked plots by latency
    latency2=latency(1:endidx);
    [~,sortedidx]=sort(latency2);
    allStacked=allStacked(:,sortedidx,:);
    event_indices(:,1:endidx)=event_indices(:,sortedidx);
    
    
    %update average
    %current_num_event=prev_num_events;
    average = (old_average*plotted_num_events + sum_windows)/(num_events);
    average_event2 = (old_average_event2*prev_num_events + sum_windows_event2)/(num_events);
    new_plot_flag=1;
    
    plotted_num_events=endidx;%Need to check if this is right
    
else
    %if no new events, plotted_num_events does not increment
    %num_events=prev_num_events;
    average = old_average;
    average_event2=old_average_event2;
    new_plot_flag=0;
end
%%


function [average,last_used_ind,allStacked,new_plot_flag, plotted_num_events,event_indices, lastSpec,lastSpec_event2,average_event2] =...
    extract_event_windows(num_events, plotted_num_events, event_indices, window_before_event,window_after_event,...
    DataAfterCAR,AmountOfData,number_electrodes ,num_freq_bands, freqs2plot,freq_band_singlestacked, ...
    old_average,old_average_event2, last_used_ind,allStacked,averages, stdevs,event_type)
%allStacked is the single stacked variable where each trial is sorted by
%latency of response. Aligned at event1, and includes event2
%old_average_event2 is updated average spectrogram aligned at event2
lastSpec=[];
lastSpec_event2=[];
window_around_event=window_before_event+window_after_event;

startidx=plotted_num_events+1;
endidx=num_events;
beginning = event_indices(1,startidx) - window_before_event;


switch event_type
    case 1
        last = event_indices(1,endidx)+ window_after_event;
        while(last>AmountOfData)&& endidx>=startidx% what if ind(num_events-1)+window/2 is still out of bounds?
            last = event_indices(1,endidx)+ window_after_event;
            endidx=num_events-1;
        end
        
        if isempty(old_average)
            old_average = zeros(number_electrodes,num_freq_bands, window_before_event*2+1);
        end
        
        sum_windows = zeros(number_electrodes,num_freq_bands, window_before_event*2+1);
        
    case 2
        last2 = event_indices(2,endidx)+ window_before_event;
        last = event_indices(1,endidx)+ window_after_event;
        while((last>AmountOfData) || (last2>AmountOfData)) && endidx>=startidx% what if ind(num_events-1)+window/2 is still out of bounds?
            last2 = event_indices(2,endidx)+ window_before_event;
            last = event_indices(1,endidx)+ window_after_event;
            endidx=num_events-1;
        end
        
        
        if isempty(old_average)
            old_average = zeros(number_electrodes,num_freq_bands, window_before_event*2+1);
            old_average_event2=zeros(number_electrodes,num_freq_bands, window_before_event*2+1);
        end
        sum_windows = zeros(number_electrodes,num_freq_bands, window_before_event*2+1);
        sum_windows_event2 = zeros(number_electrodes,num_freq_bands, window_before_event*2+1);
end


%If there are new events within the bounds of current data
%do i even need this? will endidx always be >= startidx?
if endidx>=startidx
    %loop over each event, grab corresponding window
    for i=[startidx:endidx]
        beginning = event_indices(1,i) - window_before_event;
        last = event_indices(1,i)+ window_after_event; % NOTE: need different last for stacked, see ln 761.
        %I think we should grab 2.5seconds after the first event
        %new_count=new_count+1;
        timepts = beginning:last;
        lastSpec=DataAfterCAR(:,freqs2plot, timepts);
        last_used_ind = event_indices(1,i);
        %prev_num_events=prev_num_events+1;
        lastSpec=(lastSpec-averages)./stdevs;
        sum_windows = sum_windows + lastSpec(:,:,1:window_before_event*2+1);
        
        allStacked(:,i,:)=lastSpec(:,freq_band_singlestacked,:);
        
        if event_type==2
            %Get data aligned at 2nd event
            beginning = event_indices(2,i) - window_before_event;
            last = event_indices(2,i)+ window_before_event;
            timepts = beginning:last;
            
            lastSpec_event2=(DataAfterCAR(:,freqs2plot, timepts)-averages(:,:,1:2*window_before_event+1))./stdevs(:,:,1:2*window_before_event+1);
            sum_windows_event2 = sum_windows_event2 + lastSpec_event2;
        end
    end
    
    if event_type==2
        %sort stacked plots by latency
        latency=event_indices(2,:)-event_indices(1,:);
        latency2=latency(1:endidx);
        [~,sortedidx]=sort(latency2);
        allStacked=allStacked(:,sortedidx,:);
        event_indices(:,1:endidx)=event_indices(:,sortedidx);
        average_event2 = (old_average_event2*plotted_num_events + sum_windows_event2)/(endidx);
        
    end
    %update average
    average = (old_average*plotted_num_events + sum_windows)/(endidx);
    
    new_plot_flag=1;
    plotted_num_events=endidx;%Need to check if this is right
else
    %if no new events, plotted_num_events does not increment
    %num_events=prev_num_events;
    average = old_average;
    average_event2=old_average_event2;
    new_plot_flag=0;
end


%%
function [num_events,indLastEvent,eventTimes,new_event_flag] = parseEvents(ANNewData_finalMAT, ANfinalPos, desired_ANchan,...
    desired_ANchan2, threshold, threshold2, sampling_rate, last_num_events,...
    indLastEvent, eventTimes, number_of_analog,window_before_event)
% num_events is number after rejection based off of timing of events
num_events = last_num_events;
event=(ANNewData_finalMAT(desired_ANchan,:)>threshold);
trigger1=[(diff(event)>0) 0];
detected_num_events = sum(trigger1);

if number_of_analog==1
    ind = find(trigger1,detected_num_events);
    ind(find([0 diff(ind)<(sampling_rate*1)]))=[]; %ignore events within 1sec of each other
    indLastEvent = [];
    %eventIndices(1,1:length(ind))=ind;
    eventTimes=ind/sampling_rate;
    num_events = length(ind);
else % number_of_analog==2
    
    event=(ANNewData_finalMAT(desired_ANchan2,:)>threshold2);
    trigger2=[2*(diff(event)>0) 0];
    
    triggers = trigger1+trigger2; %this takes care of when 2 events happen at the same time
    
    while indLastEvent<ANfinalPos
        % find onset of 1st events (starting from last event recorded) for the first two trials
        event1onset = find([zeros(1,indLastEvent) triggers(indLastEvent+1:ANfinalPos)]==1,2);
        if ~isempty(event1onset)
            %find 2nd event that occurs after 1st event
            event2onset =  find([zeros(1,event1onset(1)) triggers(event1onset(1)+1:ANfinalPos)]==2,1);
            
            %find previous trial's 1st and 2nd event
            latest1eventonset= find(triggers(1:event1onset(1)-1)==1);
            if ~isempty(latest1eventonset)
                latest1eventonset = latest1eventonset(end);
            else
                latest1eventonset = -1*sampling_rate;
            end
            latest2eventonset = find(triggers(1:event1onset(1))==2);
            if ~isempty(latest2eventonset)
                latest2eventonset = latest2eventonset(end);
            else
                latest2eventonset = -1*sampling_rate;
            end
            %end
        else
            event2onset=[];
        end
        
        if isempty(event2onset) || isempty(event1onset) %no events found
            break;
            % indLastEvent = ANfinalPos;
        elseif length(event1onset)~=2 % 2nd trial was not found
            break;
        elseif (event2onset-event1onset(1))<=2*sampling_rate && (event2onset-event1onset(1))>0 && ...
                (event1onset(2)-event2onset)>=1*sampling_rate && ...
                (event1onset(1)-latest2eventonset)>=1*sampling_rate &&...
                (event1onset(1)-latest1eventonset)>=1*sampling_rate &&...
                (event1onset(1)>window_before_event)% Takes care of not enough data before the first event
            % 2nd event occured within 0-2 seconds of 1st event
            % 2nd trial occured after 1 second of 2nd event
            % 1st event occured after 1 second of previous trial's 2nd event
            % 1st event occured after 1 second of previous trial's 1st event
            num_events = num_events+1;
            eventIndices(:,num_events) = [event1onset(1); event2onset; event2onset-event1onset(1)];
            indLastEvent = event1onset(1);
            
        else
            % next interation will look for events only after 2nd trial
            indLastEvent = event1onset(2)-1;
        end
    end
end

if num_events>last_num_events
    new_event_flag=1;
else
    new_event_flag=0;
end

function badChannels=detectBadChannels(data)
%badChannels=find(zscore(mean(data,2))>5);
for c=1:size(data,1)
    %p(c,:)=log10(smooth(periodogram(data(c,:))',200));
    p(c,:)=log10((periodogram(data(c,:))'));
end
badChannels=find(zscore(mean(abs((p-repmat(mean(p),256,1)))'))>3);

function inputBaselineFile_Callback(hObject, eventdata, handles)
% hObject    handle to inputBaselineFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of inputBaselineFile as text
%        str2double(get(hObject,'String')) returns contents of inputBaselineFile as a double


% --- Executes during object creation, after setting all properties.
function inputBaselineFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to inputBaselineFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in locateBaseline.
function locateBaseline_Callback(hObject, eventdata, handles)
% hObject    handle to locateBaseline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fileName,pathName] = uigetfile('baseline_stats_*.mat','Select Baseline File');
set(handles.inputBaselineFile,'string',[pathName filesep fileName]);
drawnow
guidata(hObject, handles);


function inputMRIFile_Callback(hObject, eventdata, handles)
% hObject    handle to inputMRIFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of inputMRIFile as text
%        str2double(get(hObject,'String')) returns contents of inputMRIFile as a double


% --- Executes during object creation, after setting all properties.
function inputMRIFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to inputMRIFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in locateMRI.
function locateMRI_Callback(hObject, eventdata, handles)
% hObject    handle to locateMRI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fileName,pathName] = uigetfile('.jpg','Select MRI File');
set(handles.inputMRIFile,'string',[pathName filesep fileName]);
drawnow
handles=plotMRI(handles,'initiate',[],[pathName filesep fileName])
guidata(hObject, handles);



function htkFilePath_Callback(hObject, eventdata, handles)
% hObject    handle to htkFilePath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of htkFilePath as text
%        str2double(get(hObject,'String')) returns contents of htkFilePath as a double


% --- Executes during object creation, after setting all properties.
function htkFilePath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to htkFilePath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function [ecogOut,nextRawIdx]=calcSTFT(data,sampFreq,nextRawIdx)
wl=2^nextpow2(fix(.05*sampFreq));
% for cidx=1:size(data,1)
%     [YY,tt,ff]=stft_hann_(data(cidx,:)',sampFreq,wl);
%     ecogOut.data(cidx,:,:)=(abs(YY));
% end
%[YY,tt,ff,wss,wl]=stft_multi_ECoG(data',sampFreq);
[YY,tt,ff,wss]=stft_hann_allChan(data',sampFreq,wl);
ecogOut.data=abs(YY);
ecogOut.time=tt;
ecogOut.freq=ff;
ecogOut.data=permute(ecogOut.data,[3,1,2]);
if ~isempty(wss) 
    nextRawIdx=wss(end)+wl/2+nextRawIdx-1;
else
    nextRawIdx=nextRawIdx;
end
%
% ecog.data=data;
% ecog.sampFreq=sampFreq;
% [ecog,cfs,sigma_fs]=processingHilbertTransform_filterbankGUI(ecog, sampFreq, [4 200]);
% ecog.time=(1:size(ecog.data,2))/sampFreq;
% ecogOut=ecog;
% ecogOut.data=permute(abs(ecog.data),[1 3 2]);

function handles=initiatePlots(handles,num_avgEvent_freqs2plot,num_singleEvent_freqs2plot,number_of_electrodes_total,...
    window_around_event,to_plot_grid,desired_freq_plot,number_of_analog,sampPerSecond,allStacked)

%% GET ELECTRODE NUMBERING
handles=getElectrodeNums(handles,num_avgEvent_freqs2plot,num_singleEvent_freqs2plot,number_of_electrodes_total,window_around_event,to_plot_grid);

%% CONTINUAL PLOT
handles=plotContinual(handles,'initiate');

%% AVERAGE SPECTROGRAM PLOT
to_plot=allStacked(:,1,:,:);
to_plot=squeeze(nanmean(to_plot,2));
flattened = reshape_3Ddata(to_plot, window_around_event,size(to_plot,2), to_plot_grid);
to_plot=flattened;
handles=plotAveSpec(handles,'initiate',[],[],to_plot,[],[],number_of_analog,num_avgEvent_freqs2plot);

%%
%EVENT COUNTER PLOT
handles=plotEventCounter(handles,'initiate');

%% PLOT STACKED
to_plot=mean(allStacked(:,1,desired_freq_plot,:),3);
flattened = reshape_3Ddata(to_plot, window_around_event,1, to_plot_grid);
to_plot=flattened;
handles=plotStacked(handles,'initiate',window_around_event,to_plot_grid,to_plot,number_of_electrodes_total,1)


if ~isempty(get(handles.inputMRIFile,'string'))
    handles=plotMRI(handles,'initiate',to_plot,get(handles.inputMRIFile,'string'),desired_freq_plot);
end
