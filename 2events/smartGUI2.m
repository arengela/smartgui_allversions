function varargout = smartGUI2(varargin)
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

% Last Modified by GUIDE v2.5 15-Nov-2011 16:17:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
%uicontrol('Style','pushbutton','String','collect','Position',[10 10 20 20])
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @smartGUI2_OpeningFcn, ...
    'gui_OutputFcn',  @smartGUI2_OutputFcn, ...
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
function smartGUI2_OpeningFcn(hObject, eventdata, handles, varargin)
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
function varargout = smartGUI2_OutputFcn(hObject, eventdata, handles)
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
sampling_rate = 500;
number_of_channels = 16;
ANnumber_of_channels = 4;
number_of_sec2read = 2; %******Does this need to change?
dANproj_name =  {'motormap.danin'};
sANproj_name =  {'motormap.sanin'};
integration_time = .1;
avgEvent_freqs2plot = 1:7;
singleEvent_freqs2plot = 1:7;
num_freq_bands = 7;
window_before_event=250;


% Retrieve from GUI
time2collectBaseline = str2num(get(handles.time2collectBaseline,'String'));
desired_freq_plot = str2num(get(handles.desired_freq_plot,'String'));
freq_band_singlestacked = str2num(get(handles.freq_band_singlestacked,'String'));
desired_ANchan = str2num(get(handles.desired_ANchan,'String'));
threshold = str2num(get(handles.threshold,'String'));
subj_id = get(handles.subj_id,'String');
elec_tmp = get(handles.get_total_electrodes,'String');
number_of_electrodes_total=str2num(elec_tmp{get(handles.get_total_electrodes,'Value')});
number_of_analog=get(handles.get_num_analog,'Value');
%number_of_analog=str2num(elec_tmp{get(handles.get_num_analog,'Value')});
%allStacked=zeros(number_of_electrodes_total,num_freq_bands,window_around_event+1,50);


if number_of_analog==2
    window_after_event=window_before_event+2.5*500;
else
    window_after_event=window_before_event;
end

window_around_event=window_before_event+window_after_event;
 number_of_sec2read = window_around_event/500; %******Does this need to change?
   
allStacked=zeros(number_of_electrodes_total,50,window_around_event+1);
    

% Flags
event_flag = 1;  % Average plot
calculated_baseline_flag = get(handles.calculated_baseline_flag,'Value');


if  number_of_electrodes_total == 64%get(handles.get_total_electrodes,'Value')==1
    dproj_name = {'motormap.dwav1';'motormap.dwav2';'motormap.dwav3';'motormap.dwav4'};
    sproj_name = {'motormap.swav1';'motormap.swav2';'motormap.swav3';'motormap.swav4'};
    to_plot_grid=8;
elseif number_of_electrodes_total == 256;
    dproj_name = {'motormap.dwav1';'motormap.dwav2';'motormap.dwav3';'motormap.dwav4';...
        'motormap.dwav5';'motormap.dwav6';'motormap.dwav7';'motormap.dwav8';...
        'motormap.dwav9';'motormap.dwa10';'motormap.dwa11';'motormap.dwa12';...
        'motormap.dwa13';'motormap.dwa14';'motormap.dwa15';'motormap.dwa16'};
    sproj_name = {'motormap.swav1';'motormap.swav2';'motormap.swav3';'motormap.swav4';...
        'motormap.swav5';'motormap.swav6';'motormap.swav7';'motormap.swav8';...
        'motormap.swav9';'motormap.swa10';'motormap.swa11';'motormap.swa12';...
        'motormap.swa13';'motormap.swa14';'motormap.swa15';'motormap.swa16'};
    to_plot_grid=16;
end

if number_of_analog==2
    desired_ANchan2 = str2num(get(handles.desired_ANchan2,'String'));
    threshold2 = str2num(get(handles.threshold2,'String'));
else
    desired_ANchan2 = [];
    threshold2 = [];
end



% Calculate
points_needed4baseline = sampling_rate*time2collectBaseline*number_of_channels *num_freq_bands;
number_of_points2read = sampling_rate*number_of_sec2read*number_of_channels*num_freq_bands;
ANnumber_of_points2read = sampling_rate*number_of_sec2read*ANnumber_of_channels;
num_avgEvent_freqs2plot = length(avgEvent_freqs2plot);
num_singleEvent_freqs2plot = length(singleEvent_freqs2plot);

%%
%setPlots;

%CONTINUAL PLOT
set(0,'currentfigure',handles.figure1);
set(handles.figure1,'CurrentAxes',handles.continualAverage); cla
fh2=handles.continualAverage;
title('Continual Plot')
set(fh2,'XTickLabel',[]);
set(fh2,'YTickLabel',[]);

%%
%AVERAGE SPECTROGRAM PLOT
fh=figure(2); clf% figure for average event trigger analysis
z=axes('Position',[.1 .1 .85 .85],'visible','off');
set(fh,'DefaultAxesxtick',[0.5: (window_before_event*2+1):to_plot_grid*(window_before_event*2+1)])
set(fh,'DefaultAxesytick',[0.5: num_avgEvent_freqs2plot:to_plot_grid*num_avgEvent_freqs2plot])
set(fh,'DefaultAxesyticklabel','')
set(fh,'DefaultAxesxticklabel','')
set(fh,'DefaultAxeslinewidth',1)
set(fh,'DefaultAxesgridlinestyle','-')
set(fh, 'Name','Average Event Time Frequency Plot','NumberTitle','off')
set(fh,'CurrentAxes',z);


%%
%LAST SPECTROGRAM PLOT
%fh4=figure(4); % figure for single event trigger analysis
set(0,'currentfigure',handles.figure1);
set(handles.figure1,'CurrentAxes',handles.lastSpectrogram); cla
fh4=handles.lastSpectrogram;
%set(fh4,'Position',[.1 .1 .85 .85],'visible','off');
set(fh4,'Xtick',[0.5: (window_around_event+1):to_plot_grid*(window_around_event+1)])
set(fh4,'YTick',[0.5: num_singleEvent_freqs2plot:to_plot_grid*num_singleEvent_freqs2plot])
set(fh4,'YTickLabel','')
set(fh4,'XTickLabel','')
set(fh4,'LineWidth',1)
set(fh4,'GridLineStyle','-')
title('Single Event Time Frequency Plot')
%set(fh4,'CurrentAxes',z)
%%
%EVENT COUNTER PLOT
%fh5 = figure(5);
set(0,'currentfigure',handles.figure1);
set(handles.figure1,'CurrentAxes',handles.eventCounter); cla
fh5=handles.eventCounter;
title('Number of events vs matlab loop counter')
set(fh5,'XTick',[0:100])
set(fh5,'YTick',[0:100])
set(fh5,'XTickLabel',[])
set(fh5,'YTickLabel',[])
%%
%SINGLE STACKED PLOT
fh6 = figure(6); clf %single stacked plot
z=axes('Position',[.1 .1 .85 .85],'visible','off');
set(fh6,'DefaultAxesxtick',[0.5: (window_around_event+1):to_plot_grid*(window_around_event+1)])
set(fh6,'DefaultAxeslinewidth',1)
set(fh6,'DefaultAxesyticklabel','')
set(fh6,'DefaultAxesxticklabel','')
set(fh6,'DefaultAxesgridlinestyle','-')
set(fh6, 'Name','Single Stacked','NumberTitle','off')
set(fh6,'CurrentAxes',z);

%%
%Make figure for average spectrogram of 2nd event
if number_of_analog==2
    fh_event2=figure(7);clf;      
    z=axes('Position',[.1 .1 .85 .85],'visible','off');
    set(fh_event2,'DefaultAxesxtick',[0.5: (window_before_event*2+1):to_plot_grid*(window_before_event*2+1)])
    set(fh_event2,'DefaultAxesytick',[0.5: num_avgEvent_freqs2plot:to_plot_grid*num_avgEvent_freqs2plot])
    set(fh_event2,'DefaultAxesyticklabel','')
    set(fh_event2,'DefaultAxesxticklabel','')
    set(fh_event2,'DefaultAxeslinewidth',1)
    set(fh_event2,'DefaultAxesgridlinestyle','-')
    set(fh_event2, 'Name','Average Event Time Frequency Plot','NumberTitle','off')
    set(fh_event2,'CurrentAxes',z);
end

%%
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
fig_dash_x = repmat([round((window_around_event+1)/2):window_around_event+1:to_plot_grid*+window_around_event],2,1);

if number_of_analog==2
    % coordinates for the electrode numbers
    fig_num_y2 = [1.1:num_avgEvent_freqs2plot:to_plot_grid*num_avgEvent_freqs2plot];
    fig_num_y2 = reshape(repmat(fig_num_y,[to_plot_grid 1]),number_of_electrodes_total,1);
    fig_num_x2 = [.5:(window_before_event*2+1):to_plot_grid*(window_before_event*2+1)];
    fig_num_x2 = reshape(repmat(fig_num_x,[1 to_plot_grid]),number_of_electrodes_total,1);
    fig_nums2 = {};
    for i = 1:number_of_electrodes_total
        fig_nums2 = [fig_nums2 {num2str(i)}];
    end

    % coordinates for the dash lines
    avgfig_dash_y2 = repmat([0 num_avgEvent_freqs2plot*to_plot_grid],to_plot_grid,1)';
    sigfig_dash_y2 = repmat([0 num_singleEvent_freqs2plot*to_plot_grid],to_plot_grid,1)';
    fig_dash_x2= repmat([round((window_before_event*2+1)/2):window_before_event*2+1:to_plot_grid*window_before_event*2],2,1);
end

pth = pwd;
addpath([pth filesep 'movie']);
rehash path;

%%
if calculated_baseline_flag==1
    %load baseline_stats
    load('C:\Users\Connie\Desktop\SMART\baseline_stats_EC_B.mat');
    %uiopen('load')
    medians = squeeze(medians(:,:,1));
end

bufferCounter = zeros(1,size(dproj_name,1));
ANbufferCounter = 0;
matlab_loopCounter = 0;
num_samples = 0; %used to calculate faster average
last_used_ind = 0;  %used to calculate faster average
eventRelatedAvg = [];
good_event_count = 0;
prev_num_events = 0;
indLastEvent = 0; % for 2 analogs
eventIndices = zeros(2,50); % for 2 analogs

%buffer of 5 minutes
if event_flag ==1
    DataAfterCAR=zeros(number_of_electrodes_total, num_freq_bands, 4*60*sampling_rate);
    ANNewData_finalMAT=zeros(ANnumber_of_channels, 4*60*sampling_rate);
end

DA = actxcontrol('TDevAcc.X');% Loads ActiveX methods
if DA.ConnectServer('Local')>0 %Checks to see if there is an OpenProject and Workbench load
    %     if DA.SetSysMode(2)==1%Starts OpenEx project running and acquiring data. 1 = Standby, 2 = Preview, !!MUST CHANGE TO 3 WHEN RECORDING
    if DA.GetSysMode>1 % Checks to see if the Project is running
        if not(isnan(DA.ReadTargetVEX(dproj_name{1},0,10,'F32','F64'))) %checks to see if it loaded the correct project
            for i=1:length(dproj_name)
                BufferSize(i)=DA.GetTargetSize(dproj_name{i});% returns the Buffer size (bookkeeping)
            end
            if event_flag ==1
                ANBufferSize=DA.GetTargetSize(dANproj_name{1});% returns the Analog Buffer size (bookkeeping)
                ANoldAcqPos=ANBufferSize;
            end
            oldAcqPos=BufferSize;
            
            %if DA.GetSysMode=0, then means tdt is not recording or prj not working
            while DA.GetSysMode>1
                
                tic
                %%
                %Calculate baseline
                if calculated_baseline_flag==0
                    [averages, stdevs, medians]=calculate_baseline(time2collectBaseline, ...
                        sproj_name,dproj_name, points_needed4baseline, BufferSize, DA, number_of_channels, ...
                        num_freq_bands, sampling_rate,number_of_sec2read);
                    pause(1)
                    uisave({'averages' 'stdevs' 'medians'}, ['baseline_stats_' subj_id])
                    medians = squeeze(medians(:,:,1));
                    calculated_baseline_flag=1;
                    matlab_loopCounter = 0;
                end
                %%
                %Retrieve data
                % find place in TDT buffer, check for buffer wrap
                % around, find place in MATLAB buffer
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
                
                disp([num2str(number_of_points2read/112/500) 'sec'])
                
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
                if event_flag ==1
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
                end
                
                % Read from TDT buffers and reshape to chan x freq bands x time
                NewData=readTarget(dproj_name, AcqPos, number_of_points2read, DA);
                NewDataMAT=zeros(number_of_channels*length(AcqPos),num_freq_bands,number_of_points2read/(num_freq_bands*number_of_channels));
                for i = 1:size(NewData,1)
                    ind = (i-1)*number_of_channels+1;
                    NewDataMAT( ind:(ind+number_of_channels-1),:, :) =shiftdim(reshape(reshape(NewData(i,:),...
                        num_freq_bands*number_of_channels, number_of_points2read/(num_freq_bands*number_of_channels))',...
                        number_of_points2read/(num_freq_bands*number_of_channels),number_of_channels,num_freq_bands),1);
                end
                
                %Take log(power+medians)
                tmp = repmat(medians,[1 1 (number_of_points2read/(num_freq_bands*number_of_channels))]);
                LogNewData =log(abs(NewDataMAT)+tmp+eps);
                
                %%
                %Plot continual plot
                %continual plotting
                set(0,'currentfigure',handles.figure1);
                set(handles.figure1,'CurrentAxes',fh2)
                runningAverage = squeeze(mean(LogNewData(:, desired_freq_plot, (end-integration_time*sampling_rate):end),3));
                ZscoreNewData=(runningAverage-averages(:,desired_freq_plot,1))./stdevs(:,desired_freq_plot,1);
                imagesc(real(reshape(ZscoreNewData,to_plot_grid,to_plot_grid))', [-7 7]);
                drawnow;
                %%
                %Plot event counter
                if event_flag ==1
                    % reshape ANdata by chan x time
                    ANNewDataMAT=reshape(ANNewData,ANnumber_of_channels, ANnumber_of_points2read/ANnumber_of_channels);
                    
                    % Place data into MATLAB buffers
                    % finalPos = posInNewData_AfterCAR+sampling_rate*number_of_sec2read-1;
                    finalPos = posInNewData_AfterCAR + number_of_points2read/(num_freq_bands*number_of_channels) - 1;
                    ANfinalPos = ANposInNewData + ANnumber_of_points2read/ANnumber_of_channels - 1;
                    
                    DataAfterCAR(:,:,posInNewData_AfterCAR:finalPos)=LogNewData;
                    ANNewData_finalMAT(:,ANposInNewData:ANfinalPos)=ANNewDataMAT;
                    %%
                    %find events
                    event=(ANNewData_finalMAT(desired_ANchan,:)>threshold);
                    trigger=(diff(event)>0);
                    detected_num_events = sum(trigger);
                    
                    set(0,'currentfigure',handles.figure1);
                    set(handles.figure1,'CurrentAxes',fh5);
                    hold on
                    plot(matlab_loopCounter,detected_num_events,'*');
                    grid on;
                    
                    
                    %%
                    %Plot average spectrogram
                    if detected_num_events>0 && event_flag == 1
                        
                        %%% PARSE EVENTS HERE
                        [num_events,indLastEvent,eventIndices] = ...
                            parseEvents(ANNewData_finalMAT, ANfinalPos, desired_ANchan,...
                            desired_ANchan2, threshold, threshold2, sampling_rate, prev_num_events,...
                            indLastEvent, eventIndices, number_of_analog);
                        
                        if num_events>prev_num_events
                            % num_events is #events after removing ones within 1 second, intialized to 0
                            prev_num_events = num_events;
                            newEventFlag=1;
                        else
                            newEventFlag=0;
                        end
                        
                        if  newEventFlag==1 && number_of_analog ==1
                            %%% RUN OLD AVERAGE CODE
                            %{
                            Note: Can also get rid of num_samples,
                            old_num_samples like what I did with average_event_window_2_events
                            
                            PS Note: we should test this new implementation
                            in the simulator to make sure num_samples and old_num_samples
                            is really the same as good_event_count and current_num_events
                            %}
                            [eventRelatedAvg,num_samples,last_used_ind,allStacked,new_plot_flag,good_event_count,lastSpec] = average_event_window_angela(num_events,eventIndices(1,1:num_events) , window_around_event,...
                                DataAfterCAR,finalPos(1),number_of_electrodes_total,...
                                num_avgEvent_freqs2plot, avgEvent_freqs2plot,...
                                eventRelatedAvg, num_samples, last_used_ind,allStacked,good_event_count, averages, stdevs,freq_band_singlestacked);
                            
                            
                        elseif newEventFlag==1 && number_of_analog ==2
                            
                            [eventRelatedAvg,last_used_ind,allStacked,new_plot_flag,current_num_events,lastSpec,lastSpec_event2,eventRelatedAvg2] =...
                                average_event_window_2_events(num_events, event_indices, window_around_event,window_after_event,...
                                DataAfterCAR,finalPos(1),number_electrodes ,num_freq_bands, freqs2plot,...
                                old_average, last_used_ind,allStacked,good_event_count, averages, stdevs,old_average_event2,lastSpec_event2, freq_band_singlestacked)
                        end
                        
                        if  new_plot_flag==1 & newEventFlag==1
                            
                            %Plot new event count in red
                            set(0,'currentfigure',handles.figure1);
                            set(handles.figure1,'CurrentAxes',fh5)
                            hold on
                            plot(matlab_loopCounter,good_event_count,'r*');
                            
                            %%
                            %Plot average Spectrogram(s) around event(s)
                            %always plot average spectrogram aligned
                            %around event 1
                            
                            
                            if number_of_analog ==2
                                ZscoreEventRelatedAvg=(eventRelatedAvg-averages)./stdevs;
                                to_plot=  reshape_3Ddata(ZscoreEventRelatedAvg, window_before_event*2,...
                                num_avgEvent_freqs2plot, to_plot_grid);
                                to_plot = artifactrejection(to_plot, window_before_event*2,-1.5, 1.5, .5, .8, to_plot_grid);

                                set(0,'currentfigure',fh)
                                imagesc(real(to_plot),[-3 3]);
                                text(fig_num_x2, fig_num_y2, fig_nums)
                                line(fig_dash_x2,avgfig_dash_y2,'LineStyle','--','linewidth',1','color',[0 0 0]);
                                grid on;
                                drawnow;
                                
                                ZscoreEventRelatedAvg=(eventRelatedAvg2-averages)./stdevs;
                                to_plot=  reshape_3Ddata(ZscoreEventRelatedAvg, window_before_event*2,...
                                    num_avgEvent_freqs2plot, to_plot_grid);
                                to_plot = artifactrejection(to_plot, window_before_event*2,-1.5, 1.5, .5, .8, to_plot_grid);
                                
                                %Plot average spectrogram aligned at
                                %event 2
                                set(0,'currentfigure',fh_event2)
                                imagesc(real(to_plot),[-3 3]);
                                text(fig_num_x2, fig_num_y2, fig_nums)
                                line(fig_dash_x2,avgfig_dash_y2,'LineStyle','--','linewidth',1','color',[0 0 0]);
                                grid on;
                                drawnow;
                            else
                                ZscoreEventRelatedAvg=(eventRelatedAvg-averages)./stdevs;
                                
                                to_plot=  reshape_3Ddata(ZscoreEventRelatedAvg, window_around_event,...
                                    num_avgEvent_freqs2plot, to_plot_grid);
                                to_plot = artifactrejection(to_plot, window_around_event,-1.5, 1.5, .5, .8, to_plot_grid);
                                
                                set(0,'currentfigure',fh)
                                imagesc(real(to_plot),[-3 3]);
                                text(fig_num_x, fig_num_y, fig_nums)
                                line(fig_dash_x,avgfig_dash_y,'LineStyle','--','linewidth',1','color',[0 0 0]);
                                grid on;
                                drawnow;
                            end                           
                            
                            
                            %%
                            %Plot last spectrogram
                            ZscoreLastEvent=lastSpec;
                            
                            to_plot = reshape_3Ddata(ZscoreLastEvent, window_around_event,...
                                num_singleEvent_freqs2plot, to_plot_grid);
                            
                            set(0,'currentfigure',handles.figure1);
                            set(handles.figure1,'CurrentAxes',fh4)
                            imagesc(real(to_plot),[-7 7]);
                            text(fig_num_x, fig_num_y, fig_nums)
                            line(fig_dash_x,sigfig_dash_y,'LineStyle','--','linewidth',1','color',[0 0 0]);
                            set(fh4,'Xtick',[0.5: (window_around_event+1):to_plot_grid*(window_around_event+1)])
                            set(fh4,'YTick',[0.5: num_singleEvent_freqs2plot:to_plot_grid*num_singleEvent_freqs2plot])
                            set(fh4,'YTickLabel','')
                            set(fh4,'XTickLabel','')
                            set(fh4,'LineWidth',1)
                            set(fh4,'GridLineStyle','-')
                            grid on;
                            drawnow;
                            %%
                            %plot single stacked
                            
                            to_plot=allStacked(:,1:good_event_count,:);
                            
                            set(fh6,'DefaultAxesytick',[0.5: good_event_count:to_plot_grid*good_event_count])
                            set(0,'currentfigure',fh6);
                            
                            %Will reshape functions work with longer
                            %window after first event?
                            flattened = reshape_3Ddata(to_plot, window_around_event,good_event_count, to_plot_grid);                            
                            
                            to_plot = artifactrejection(flattened, window_around_event,-5.5, 5.5, .5, .8, to_plot_grid);
                            imagesc(real(to_plot),[-7 7]);
                            if good_event_count >1
                                fig_num_y_stacked = [1.1:good_event_count:to_plot_grid*good_event_count];
                                fig_num_y_stacked = reshape(repmat(fig_num_y_stacked, [to_plot_grid 1]), number_of_electrodes_total, 1);
                                text(fig_num_x, fig_num_y_stacked, fig_nums)
                                stackedfig_dash_y = repmat([0 good_event_count*to_plot_grid],to_plot_grid,1)';
                                line(fig_dash_x,stackedfig_dash_y,'LineStyle','--','linewidth',1','color',[0 0 0]);
                            end
                            grid on
                            
                        end
                        %end
                    end
                    
                    ANoldAcqPos=ANAcqPos;
                end
                
                matlab_loopCounter=matlab_loopCounter+1;
                oldAcqPos=AcqPos;
                
                if toc< integration_time*sampling_rate
                    pause(.1)
                end
            end
            
            
        else
            msgbox('Incorrect OpenEx Project')
        end
    else
        msgbox('OpenEx project Failed To Run')
    end
    %     end
else
    msgbox('OpenEx project not loaded reload OpenEx project and restart MATLAB script')
end
DA.CloseConnection
disp('Saving figures...')
saveas(fh,[subj_id '_avg.jpg'],'jpg')
saveas(fh6,[subj_id '_stacked.jpg'],'jpg')

play_movie = input('Play movie (1/0): ');
while play_movie
    save_movie = 0;
    makeMovie(subj_id, ZscoreEventRelatedAvg,to_plot_grid, number_of_electrodes_total, save_movie);
    play_movie = input('Play movie (1/0): ');
end

uisave({'ZscoreEventRelatedAvg' }, ['ZscoreEventRelatedAvg' subj_id])
% disp('Save data?')
% uisave({'DataAfterCAR' 'ANNewData_finalMAT'}, ['data_' subj_id])

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

%%%%%%
% function NewData_AfterCAR = commonAveraging(NewDataMAT, good_channels4CAR, number_of_channels)
% %common averaging of EEG
% CAR=mean(NewDataMAT(good_channels4CAR,:),1);
% CAR_Mat=repmat(CAR,number_of_channels,1);
% NewDataMAT_AfterCAR=NewDataMAT-CAR_Mat;
%
% %increase precision, decrease bit
% NewData_AfterCAR=1000*NewDataMAT_AfterCAR;

%%%%%%

function [average,num_samples,last_used_ind,allStacked,new_plot_flag,current_num_events,lastSpec] = average_event_window_angela(num_events, ind, window_around_event,...
    DataAfterCAR,AmountOfData,number_electrodes ,num_freq_bands, freqs2plot,...
    old_average, old_num_samples, last_used_ind,allStacked,good_event_count, averages, stdevs,freq_band_singlestacked)
lastSpec=[];
%allocate memory, set up counter
current_num_events=num_events;
% average = zeros(number_electrodes,num_freq_bands, window_around_event+1);
sum_windows = zeros(number_electrodes,num_freq_bands, window_around_event+1);
num_samples = 0;

if isempty(old_average)
    old_average = zeros(number_electrodes,num_freq_bands, window_around_event+1);
end
if isempty(find(ind==last_used_ind))
    x=0;
else
    x=find(ind==last_used_ind);
end
%loop over each event, grab corresponding window
new_count=0;
for i=(x+1):current_num_events
    beginning = ind(i) - window_around_event/2;
    last = ind(i)+ window_around_event/2;
    if beginning <1 || last>AmountOfData
        new_count=0;
        continue; %since can't add different sized window, just ignore
    end
    new_count=1;
    timepts = beginning:last;
    lastSpec=DataAfterCAR(:,freqs2plot, timepts);
    sum_windows = sum_windows + DataAfterCAR(:,freqs2plot, timepts);
    num_samples = num_samples + 1;
    last_used_ind = ind(i);
    current_num_events=current_num_events+new_count;
   lastSpec=(lastSpec-averages)./stdevs;
   allStacked(:,current_num_events,:)=lastSpec(:,freq_band_singlestacked,:);
    new_plot_flag=1;
end;
if num_samples>0
    num_samples = num_samples+old_num_samples;
    average = (old_average*old_num_samples + sum_windows)/(num_samples);
else %if no new events
    average = old_average;
    num_samples = old_num_samples;
    new_plot_flag=0;
end
%%
function [average,last_used_ind,allStacked,new_plot_flag,current_num_events,lastSpec,lastSpec_event2,average_event2] =...
    average_event_window_2_events(num_events, event_indices, window_around_event,window_after_event,...
    DataAfterCAR,AmountOfData,number_electrodes ,num_freq_bands, freqs2plot,...
    old_average, last_used_ind,allStacked,good_event_count, averages, stdevs,old_average_event2,lastSpec_event2, freq_band_singlestacked)
%allStacked is the single stacked variable where each trial is sorted by
%latency of response. Aligned at event1, and includes event2
%old_average_event2 is updated average spectrogram aligned at event2

window_before_event=window_around_event/2;

ind=event_indices(1,:);
latency=event_indices(2,:)-ind;
%allocate memory, set up counter
current_num_events=good_event_count;
sum_windows = zeros(number_electrodes,num_freq_bands, window_before_event*2+1);
sum_windows_event2 = zeros(number_electrodes,num_freq_bands, window_before_event*2+1);

if isempty(old_average)
    old_average = zeros(number_electrodes,num_freq_bands, window_before_event*2+1);
    old_average_event2=zeros(number_electrodes,num_freq_bands, window_before_event*2+1);
end 

%****why use find? aren't new indices just added to the end of the index
%list?
if isempty(find(ind==last_used_ind))==1
    x=0;
else
    x=find(ind==last_used_ind);
end

%Finds event segments within bounds of available data
startidx=x+1;
beginning = ind(startidx) - window_before_event;
while(beginning<1) & startidx<length(ind) %what if ind(x+2) -window/2 is still out of bounds?
    startidx=x+1;
    beginning = ind(startidx) - window_before_event;
end

endidx=length(ind);
last = event_indices(2,endidx)+ window_before_event;
last2 = ind(endidx)+ window_after_event;
while((last>AmountOfData) & (last2>AmountOfData)) & endidx>=startidx % what if ind(num_events-1)+window/2 is still out of bounds?
    endidx=num_events-1;
    last = event_indices(2,endidx)+ window_before_event;
    last2 = ind(endidx)+ window_after_event;
end

%If there are new events within the bounds of current data
%do i even need this? will endidx always be >= startidx?
if endidx>=startidx    
    %loop over each event, grab corresponding window
    new_count=0;
    for i=[startidx:endidx]
        beginning = ind(i) - window_before_event;
        last = ind(i)+ window_after_event; % NOTE: need different last for stacked, see ln 761.
        %I think we should grab 2.5seconds after the first event
        new_count=new_count+1;
        timepts = beginning:last;
        lastSpec=DataAfterCAR(:,freqs2plot, timepts);
        sum_windows = sum_windows + lastSpec(:,:,1:window_before_event*2);
        last_used_ind = ind(i);
        current_num_events=good_event_count+new_count;
        lastSpec=(lastSpec-averages)./stdevs;
        allStacked(:,current_num_events,:)=lastSpec(:,freq_band_singlestacked,:);
        
        %Get data aligned at 2nd event
        beginning = event_indices(2,i) - window_before_event;
        last = event_indices(2,i)+ window_before_event;
        timepts = beginning:last;
        
        lastSpec_event2=(DataAfterCAR(:,freqs2plot, timepts)-averages)./stdevs;
        sum_windows_event2 = sum_windows_event2 + lastSpec_event2;
        
    end;
    
    %sort stacked plots by latency
    latency2=latency(1:endidx);
    [~,sortedidx]=sort(latency2);
    allStacked=allStacked(:,sortedidx,:);
    
end

%Update averages
if current_num_events~=good_event_count
    average = (old_average*good_event_count + sum_windows)/(current_num_events);
    average_event2 = (old_average_event2*good_event_count + sum_windows_event2)/(current_num_events);    
    new_plot_flag=1;
else %if no new events
    average = old_average;
    average_event2=old_average_event2;
    new_plot_flag=0;
end

%%
function [num_events,indLastEvent,eventIndices] = parseEvents(ANNewData_finalMAT, ANfinalPos, desired_ANchan,...
    desired_ANchan2, threshold, threshold2, sampling_rate, prev_num_events,...
    indLastEvent, eventIndices, number_of_analog)
% num_events is number after rejection based off of timing of events
num_events = prev_num_events;
event=(ANNewData_finalMAT(desired_ANchan,:)>threshold);
trigger1=[(diff(event)>0) 0];
detected_num_events = sum(trigger1);

if number_of_analog==1
    ind = find(trigger1,detected_num_events);
    ind(find([0 diff(ind)<(sampling_rate*1)]))=[]; %ignore events within 1sec of each other
    indLastEvent = [];
    eventIndices(1,1:length(ind))=ind;
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
        end
        
        if isempty(event2onset) || isempty(event1onset) %no events found
            indLastEvent = ANfinalPos;
        elseif length(event1onset)~=2 % 2nd trial was not found
            break;
        elseif (event2onset-event1onset(1))<=2*sampling_rate && (event2onset-event1onset(1))>0 && ...
                (event1onset(2)-event2onset)>=1*sampling_rate && ...
                (event1onset(1)-latest2eventonset)>=1*sampling_rate &&...
                (event1onset(1)-latest1eventonset)>=1*sampling_rate
            % 2nd event occured within 0-2 seconds of 1st event
            % 2nd trial occured after 1 second of 2nd event
            % 1st event occured after 1 second of previous trial's 2nd event
            % 1st event occured after 1 second of previous trial's 1st event
            num_events = num_events+1;
            eventIndices(:,num_events) = [event1onset(1); event2onset];
            indLastEvent = event1onset(1);
            
        else
            % next interation will look for events only after 2nd trial
            indLastEvent = event1onset(2)-1;
        end
    end
end


%%
%{
function [average,num_samples,last_used_ind,allStacked,new_plot_flag,current_num_events] = average_event_window_2_events(num_events, ind, window_around_event,...
    DataAfterCAR,AmountOfData,number_electrodes ,num_freq_bands, freqs2plot,...
    old_average, old_num_samples, last_used_ind,allStacked,prev_num_events, averages, stdevs)

%allocate memory, set up counter
current_num_events=prev_num_events;
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
%loop over each event, grab corresponding window
new_count=0;
for i=(x+1):num_events
    beginning = ind(i) - window_around_event/2;
    last = ind(i)+ window_around_event/2;
    if beginning <1 || last>AmountOfData
        new_count=0;
        continue; %since can't add different sized window, just ignore
    end
    new_count=1;
    timepts = beginning:last;
    lastSpec=DataAfterCAR(:,freqs2plot, timepts);
    sum_windows = sum_windows + DataAfterCAR(:,freqs2plot, timepts);
    num_samples = num_samples + 1;
    last_used_ind = ind(i);
    current_num_events=current_num_events+new_count;
    allStacked(:,:,:,current_num_events)=(lastSpec-averages)./stdevs;
    new_plot_flag=1;
end;
if num_samples>0
    num_samples = num_samples+old_num_samples;
    average = (old_average*old_num_samples + sum_windows)/(num_samples);
else %if no new events
    average = old_average;
    num_samples = old_num_samples;
    new_plot_flag=0;
end
%}

%%


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



%%
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




%%
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

%%
function [averages, stdevs, medians]=calculate_baseline(time2collectBaseline, ...
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
for i = 1:size(BaselineData,1)
    BaselineDataMAT= [BaselineDataMAT; shiftdim(reshape(reshape(BaselineData(i,:),...
        num_freq_bands*num_channels, sampling_rate*time2collectBaseline)',...
        sampling_rate*time2collectBaseline,num_channels,num_freq_bands),1)];
end

medians = median(BaselineDataMAT,3);
logBaselineDataMAT = log(BaselineDataMAT+repmat(medians,[1 1 size(BaselineDataMAT,3)])+eps);

medians=repmat(medians,[1 1 (sampling_rate*number_of_sec2read)]);
% averages = repmat(mean(logBaselineDataMAT,3),[1 1 (sampling_rate*number_of_sec2read+1)]);
% stdevs = repmat(std(logBaselineDataMAT,0,3),[1 1 (sampling_rate*number_of_sec2read+1)]);
averages = repmat(mean(logBaselineDataMAT,3),[1 1 (sampling_rate*1+1)]);
stdevs = repmat(std(logBaselineDataMAT,0,3),[1 1 (sampling_rate*1+1)]);
waitbar(1,h,'Baseline Statistics Gathered!')
pause(1)
close(h)

%%

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

%%


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
%%
function Simulator_Callback(hObject, eventdata, handles)

a = pwd;

f=fields(handles);
for i=[24:32 34:length(fields(handles))]
    eval(sprintf('%s=deal(handles.(f{i}));',f{i}));
end

desired_ANchan = 1;
desired_ANchan2 = 2; % ADD TO GUI 


%[filename, pth] = uigetfile('*', 'Load start index and threshold');
%load([pth filename])
gdat = Envelope(:,:,start_ind:end);
analog_dat = AnalogData(:,start_ind:end);

% gdat = Envelope;
% analog_dat = AnalogData;
clear Envelope AnalogData RawData
%%


% Set variables
sampling_rate = 500;
number_of_channels = 16;
ANnumber_of_channels = 4;
number_of_sec2read = 2;
dANproj_name =  {'motormap.danin'};
sANproj_name =  {'motormap.sanin'};
integration_time = .1;
avgEvent_freqs2plot = 1:7;
singleEvent_freqs2plot = 1:7;
num_freq_bands = 7;
window_around_event = 500;
detected_num_events=0;
% Flags
event_flag = 1;  % Average plot
single_event_flag = 1;
single_stacked_flag = 1;
calculated_baseline_flag = 1;%get(handles.calculated_baseline_flag,'Value');
newEventFlag=0;
num_events=0;

% Retrieve from GUI
time2collectBaseline = str2num(get(handles.time2collectBaseline,'String'));
desired_freq_plot = str2num(get(handles.desired_freq_plot,'String'));
freq_band_singlestacked = str2num(get(handles.freq_band_singlestacked,'String'));
desired_ANchan = str2num(get(handles.desired_ANchan,'String'));
threshold = handles.threshold;
subj_id = get(handles.subj_id,'String');
number_of_electrodes_total=256;
%stackedEventsAll=zeros(number_of_electrodes_total,50,window_around_event+1);
%allStacked=zeros(number_of_electrodes_total,num_freq_bands,window_around_event+1,50);
allStacked=zeros(number_of_electrodes_total,50,window_around_event+1);
desired_ANchan = 1;
desired_ANchan2 = 2; % ADD TO GUI 
num_events = 1;
indLastEvent = 0;
eventIndices = zeros(2,50); % 50 = max number of events expected
numeventFLAG = 2; % 1 or 2  ADD TO GUI 

if  number_of_electrodes_total == 64%get(handles.get_total_electrodes,'Value')==1
    dproj_name = {'motormap.dwav1';'motormap.dwav2';'motormap.dwav3';'motormap.dwav4'};
    sproj_name = {'motormap.swav1';'motormap.swav2';'motormap.swav3';'motormap.swav4'};
    to_plot_grid=8;
elseif number_of_electrodes_total == 256;
    dproj_name = {'motormap.dwav1';'motormap.dwav2';'motormap.dwav3';'motormap.dwav4';...
        'motormap.dwav5';'motormap.dwav6';'motormap.dwav7';'motormap.dwav8';...
        'motormap.dwav9';'motormap.dwa10';'motormap.dwa11';'motormap.dwa12';...
        'motormap.dwa13';'motormap.dwa14';'motormap.dwa15';'motormap.dwa16'};
    sproj_name = {'motormap.swav1';'motormap.swav2';'motormap.swav3';'motormap.swav4';...
        'motormap.swav5';'motormap.swav6';'motormap.swav7';'motormap.swav8';...
        'motormap.swav9';'motormap.swa10';'motormap.swa11';'motormap.swa12';...
        'motormap.swa13';'motormap.swa14';'motormap.swa15';'motormap.swa16'};
    to_plot_grid=16;
end

% Calculate
points_needed4baseline = sampling_rate*time2collectBaseline*number_of_channels *num_freq_bands;
number_of_points2read = sampling_rate*number_of_sec2read*number_of_channels*num_freq_bands;
ANnumber_of_points2read = sampling_rate*number_of_sec2read*ANnumber_of_channels;
num_avgEvent_freqs2plot = length(avgEvent_freqs2plot);
num_singleEvent_freqs2plot = length(singleEvent_freqs2plot);
%%
fh=figure(2);
fh4=handles.lastSpectrogram;
fh5=handles.eventCounter;
fh2=handles.continualAverage;
fh6 = figure(6);
%setPlots;

%%
%CONTINUAL PLOT
set(0,'currentfigure',handles.figure1);
set(handles.figure1,'CurrentAxes',handles.continualAverage); cla
fh2=handles.continualAverage;
title('Continual Plot')
set(fh2,'XTickLabel',[]);
set(fh2,'YTickLabel',[]);

%%
%AVERAGE SPECTROGRAM PLOT
fh=figure(2); clf% figure for average event trigger analysis
z=axes('Position',[.1 .1 .85 .85],'visible','off');
set(fh,'DefaultAxesxtick',[0.5: (window_around_event+1):to_plot_grid*(window_around_event+1)])
set(fh,'DefaultAxesytick',[0.5: num_avgEvent_freqs2plot:to_plot_grid*num_avgEvent_freqs2plot])
set(fh,'DefaultAxesyticklabel','')
set(fh,'DefaultAxesxticklabel','')
set(fh,'DefaultAxeslinewidth',1)
set(fh,'DefaultAxesgridlinestyle','-')
set(fh, 'Name','Average Event Time Frequency Plot','NumberTitle','off')
set(fh,'CurrentAxes',z);


%%
%LAST SPECTROGRAM PLOT
%fh4=figure(4); % figure for single event trigger analysis
set(0,'currentfigure',handles.figure1);
set(handles.figure1,'CurrentAxes',handles.lastSpectrogram); cla
fh4=handles.lastSpectrogram;
%set(fh4,'Position',[.1 .1 .85 .85],'visible','off');
set(fh4,'Xtick',[0.5: (window_around_event+1):to_plot_grid*(window_around_event+1)])
set(fh4,'YTick',[0.5: num_singleEvent_freqs2plot:to_plot_grid*num_singleEvent_freqs2plot])
set(fh4,'YTickLabel','')
set(fh4,'XTickLabel','')
set(fh4,'LineWidth',1)
set(fh4,'GridLineStyle','-')
title('Single Event Time Frequency Plot')
%set(fh4,'CurrentAxes',z)
%%
%EVENT COUNTER PLOT
%fh5 = figure(5);
set(0,'currentfigure',handles.figure1);
set(handles.figure1,'CurrentAxes',handles.eventCounter); cla
fh5=handles.eventCounter;
title('Number of events vs matlab loop counter')
set(fh5,'XTick',[0:100])
set(fh5,'YTick',[0:100])
set(fh5,'XTickLabel',[])
set(fh5,'YTickLabel',[])
%%
%SINGLE STACKED PLOT
fh6 = figure(6); clf %single stacked plot
z=axes('Position',[.1 .1 .85 .85],'visible','off');
set(fh6,'DefaultAxesxtick',[0.5: (window_around_event+1):to_plot_grid*(window_around_event+1)])
set(fh6,'DefaultAxeslinewidth',1)
set(fh6,'DefaultAxesyticklabel','')
set(fh6,'DefaultAxesxticklabel','')
set(fh6,'DefaultAxesgridlinestyle','-')
set(fh6, 'Name','Single Stacked','NumberTitle','off')
set(fh6,'CurrentAxes',z);

%%
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
fig_dash_x = repmat([round((window_around_event+1)/2):window_around_event+1:to_plot_grid*window_around_event],2,1);

%%

pth = pwd;
%addpath([pth filesep 'movie']);
%rehash path;


% Variables
calculated_baseline_flag = 1;
desired_ANchan=1;

event_flag =1;
single_event_flag=1;
single_stacked_flag=1;

sampling_rate = 500;
time2collectBaseline = 30;
number_of_sec2read = 1;
window_around_event = sampling_rate;
num_singleEvent_freqs2plot = 7;
singleEvent_freqs2plot = 1:7;
num_avgEvent_freqs2plot = 7;
avgEvent_freqs2plot= 1:7;
desired_freq_plot = 6; % normally set to 5
integration_time = .1;
number_of_electrodes_total = size(gdat,1);
freq_band_singlestacked = 6;
to_plot_grid = sqrt(size(gdat,1));  % if 256, will plot 16x16 grid

%%
points_needed4baseline = sampling_rate*time2collectBaseline;
number_of_points2read = sampling_rate*number_of_sec2read;

index = 1:number_of_points2read:length(gdat);


fig_num_x = [10:(window_around_event+1):to_plot_grid*(window_around_event+1)];
fig_num_x = reshape(repmat(fig_num_x,[1 to_plot_grid]),number_of_electrodes_total,1);
fig_nums = {};
for i = 1:number_of_electrodes_total
    fig_nums = [fig_nums {num2str(i)}];
end

fig_num_x_cont = [0:to_plot_grid-1]+.55;
fig_num_x_cont = reshape(repmat(fig_num_x_cont,[1 to_plot_grid]),number_of_electrodes_total,1);
fig_num_y_cont = [0:to_plot_grid-1]+.63;
fig_num_y_cont = reshape(repmat(fig_num_y_cont,[to_plot_grid 1]),number_of_electrodes_total,1);


fig_dash_x = repmat([round((window_around_event+1)/2):window_around_event+1:to_plot_grid*window_around_event],2,1);




if calculated_baseline_flag==0
    figure(6);
    subplot(2,1,1); plot(squeeze(gdat(1,1,:)));
    subplot(2,1,2); plot(analog_dat(desired_ANchan,:));
    disp(['Points needed for 1 sec of data: ' num2str(sampling_rate)])
    start = input('Start of baseline index: ');
    final = input('End of baseline index: ');
    BaselineDataMAT = gdat(:,:,start:final);
    [averages, stdevs, medians]=calculate_baseline(BaselineDataMAT, ...
        sampling_rate, number_of_sec2read);
end

m = medians(:,:,1);
m = repmat(m, [1 1 length(gdat)]);
logGdat = log(abs(gdat)+m+eps);

matlab_loopCounter = 0;
eventRelatedAvg=[];
num_samples=0;
last_used_ind=0;%%
profile on;
good_event_count=0;
desired_ANchan = 1;
desired_ANchan2 = 2; % ADD TO GUI 
threshold = .5;
threshold2 = .7; % ADD TO GUI 
sampling_rate = 3;
num_events = 1;
indLastEvent = 0;
eventIndices = zeros(2,50); % 50 = max number of events expected
numeventFLAG = 2; % 1 or 2  ADD TO GUI 
prev_num_events=0;
number_of_analog=2;
new_plot_flag=0;
sampling_rate=500;
for i = 1:length(index)-1
    LogNewData = logGdat(:,:,index(i):index(i+1)-1);
    
    %continual plotting
    
    %%
    %Plot continual plot
    %continual plotting
    set(0,'currentfigure',handles.figure1);
    set(handles.figure1,'CurrentAxes',fh2)
    runningAverage = squeeze(mean(LogNewData(:, desired_freq_plot, (end-integration_time*sampling_rate):end),3));
    ZscoreNewData=(runningAverage-averages(:,desired_freq_plot,1))./stdevs(:,desired_freq_plot,1);
    imagesc(real(reshape(ZscoreNewData,to_plot_grid,to_plot_grid))', [-7 7]);
    drawnow;
    %%
    %Plot event counter
    if event_flag ==1
        % reshape ANdata by chan x time
        ANNewData_finalMAT= analog_dat(:,1:index(i+1)-1);
         DataAfterCAR = logGdat(:,:,1:index(i+1)-1);
        ANfinalPos=index(i+1)-1;
        % Place data into MATLAB buffers
        % finalPos = posInNewData_AfterCAR+sampling_rate*number_of_sec2read-1;
        %finalPos = posInNewData_AfterCAR + number_of_points2read/(num_freq_bands*number_of_channels) - 1;
        %ANfinalPos = ANposInNewData + ANnumber_of_points2read/ANnumber_of_channels - 1;
        
        %DataAfterCAR(:,:,posInNewData_AfterCAR:finalPos)=LogNewData;
        %ANNewData_finalMAT(:,ANposInNewData:ANfinalPos)=ANNewDataMAT;
        %%
        %find events
        prev_num_events=detected_num_events;
        event=(ANNewData_finalMAT(desired_ANchan,:)>threshold);
        trigger=(diff(event)>0);
        detected_num_events = sum(trigger);
        eventIndices(1,1:detected_num_events)=event;
        
        set(0,'currentfigure',handles.figure1);
        set(handles.figure1,'CurrentAxes',fh5);
        hold on
        plot(matlab_loopCounter,detected_num_events,'*');
        grid on;        

        
        %%
        %Plot average spectrogram
        if detected_num_events>0 && event_flag == 1
            
            %%% PARSE EVENTS HERE
            [num_events,indLastEvent,eventIndices] = ...
                parseEvents(ANNewData_finalMAT, ANfinalPos, desired_ANchan,...
                desired_ANchan2, threshold, threshold2, sampling_rate, prev_num_events,...
                indLastEvent, eventIndices, number_of_analog);
            
            if num_events>prev_num_events
                % num_events is #events after removing ones within 1 second, intialized to 0
                prev_num_events = num_events;
                newEventFlag=1;
            else
                newEventFlag=0;
            end
            
            if  newEventFlag==1 && number_of_analog ==1
                %%% RUN OLD AVERAGE CODE
                %{
                            Note: Can also get rid of num_samples,
                            old_num_samples like what I did with average_event_window_2_events
                            
                            PS Note: we should test this new implementation
                            in the simulator to make sure num_samples and old_num_samples
                            is really the same as good_event_count and current_num_events
                %}
                [eventRelatedAvg,num_samples,last_used_ind,allStacked,new_plot_flag,good_event_count] = average_event_window_angela(num_events, ind, window_around_event,...
                    DataAfterCAR,finalPos(1),number_of_electrodes_total,...
                    num_avgEvent_freqs2plot, avgEvent_freqs2plot,...
                    eventRelatedAvg, num_samples, last_used_ind,allStacked,good_event_count, averages, stdevs,freq_band_singlestacked);
                
                
            elseif newEventFlag==1 && number_of_analog ==2
                
                
                [eventRelatedAvg,last_used_ind,allStacked,new_plot_flag,current_num_events,lastSpec,lastSpec_event2,eventRelatedAvg2] =...
                    average_event_window_2_events(num_events, event_indices, window_around_event,window_after_event,...
                    DataAfterCAR,finalPos(1),number_electrodes ,num_freq_bands, freqs2plot,...
                    old_average, last_used_ind,allStacked,good_event_count, averages, stdevs,old_average_event2,lastSpec_event2, freq_band_singlestacked)
            end
            
            if  new_plot_flag==1 & newEventFlag==1
                
                %Plot new event count in red
                set(0,'currentfigure',handles.figure1);
                set(handles.figure1,'CurrentAxes',fh5)
                hold on
                plot(matlab_loopCounter,good_event_count,'r*');
                
                %%
                %Plot average Spectrogram(s) around event(s)
                %always plot average spectrogram aligned
                %around event 1
                
                
                if number_of_analog ==2
                    ZscoreEventRelatedAvg=(eventRelatedAvg-averages)./stdevs;
                    to_plot=  reshape_3Ddata(ZscoreEventRelatedAvg, window_before_event*2,...
                        num_avgEvent_freqs2plot, to_plot_grid);
                    to_plot = artifactrejection(to_plot, window_before_event*2,-1.5, 1.5, .5, .8, to_plot_grid);
                    
                    set(0,'currentfigure',fh)
                    imagesc(real(to_plot),[-3 3]);
                    text(fig_num_x2, fig_num_y2, fig_nums)
                    line(fig_dash_x2,avgfig_dash_y2,'LineStyle','--','linewidth',1','color',[0 0 0]);
                    grid on;
                    drawnow;
                    
                    ZscoreEventRelatedAvg=(eventRelatedAvg2-averages)./stdevs;
                    to_plot=  reshape_3Ddata(ZscoreEventRelatedAvg, window_before_event*2,...
                        num_avgEvent_freqs2plot, to_plot_grid);
                    to_plot = artifactrejection(to_plot, window_before_event*2,-1.5, 1.5, .5, .8, to_plot_grid);
                    
                    %Plot average spectrogram aligned at
                    %event 2
                    set(0,'currentfigure',fh_event2)
                    imagesc(real(to_plot),[-3 3]);
                    text(fig_num_x2, fig_num_y2, fig_nums)
                    line(fig_dash_x2,avgfig_dash_y2,'LineStyle','--','linewidth',1','color',[0 0 0]);
                    grid on;
                    drawnow;
                else
                    ZscoreEventRelatedAvg=(eventRelatedAvg-averages)./stdevs;
                    
                    to_plot=  reshape_3Ddata(ZscoreEventRelatedAvg, window_around_event,...
                        num_avgEvent_freqs2plot, to_plot_grid);
                    to_plot = artifactrejection(to_plot, window_around_event,-1.5, 1.5, .5, .8, to_plot_grid);
                    
                    set(0,'currentfigure',fh)
                    imagesc(real(to_plot),[-3 3]);
                    text(fig_num_x, fig_num_y, fig_nums)
                    line(fig_dash_x,avgfig_dash_y,'LineStyle','--','linewidth',1','color',[0 0 0]);
                    grid on;
                    drawnow;
                end
                
                
                %%
                %Plot last spectrogram
                ZscoreLastEvent=tmp;
                
                to_plot = reshape_3Ddata(ZscoreLastEvent, window_around_event,...
                    num_singleEvent_freqs2plot, to_plot_grid);
                
                set(0,'currentfigure',handles.figure1);
                set(handles.figure1,'CurrentAxes',fh4)
                imagesc(real(to_plot),[-7 7]);
                text(fig_num_x, fig_num_y, fig_nums)
                line(fig_dash_x,sigfig_dash_y,'LineStyle','--','linewidth',1','color',[0 0 0]);
                set(fh4,'Xtick',[0.5: (window_around_event+1):to_plot_grid*(window_around_event+1)])
                set(fh4,'YTick',[0.5: num_singleEvent_freqs2plot:to_plot_grid*num_singleEvent_freqs2plot])
                set(fh4,'YTickLabel','')
                set(fh4,'XTickLabel','')
                set(fh4,'LineWidth',1)
                set(fh4,'GridLineStyle','-')
                grid on;
                drawnow;
                %%
                %plot single stacked
                
                to_plot=allStacked(:,1:good_event_count,:);
                
                set(fh6,'DefaultAxesytick',[0.5: good_event_count:to_plot_grid*good_event_count])
                set(0,'currentfigure',fh6);
                
                %Will reshape functions work with longer
                %window after first event?
                flattened = reshape_3Ddata(to_plot, window_around_event,good_event_count, to_plot_grid);
                
                to_plot = artifactrejection(flattened, window_around_event,-5.5, 5.5, .5, .8, to_plot_grid);
                imagesc(real(to_plot),[-7 7]);
                if good_event_count >1
                    fig_num_y_stacked = [1.1:good_event_count:to_plot_grid*good_event_count];
                    fig_num_y_stacked = reshape(repmat(fig_num_y_stacked, [to_plot_grid 1]), number_of_electrodes_total, 1);
                    text(fig_num_x, fig_num_y_stacked, fig_nums)
                    stackedfig_dash_y = repmat([0 good_event_count*to_plot_grid],to_plot_grid,1)';
                    line(fig_dash_x,stackedfig_dash_y,'LineStyle','--','linewidth',1','color',[0 0 0]);
                end
                grid on
                
            end
            %end
        end
        
    end
    
end
  
    matlab_loopCounter = matlab_loopCounter+1;
    pause(.02)
profile off
keyboard
%uisave({'avg_to_plot', 'finalstackeddata', 'sing_to_plot'}, 'to_plot')

%simulator_motormapping_SMART_modified(handles)
%{
DataAfterCAR,AmountOfData,number_electrodes ,num_freq_bands, freqs2plot,...
    old_average, old_num_samples, last_used_ind,allStacked,prev_num_events, averages, stdevs

%allocate memory, set up counter
tmp=[];
current_num_events=prev_num_events;
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
%loop over each event, grab corresponding window
new_count=0;
for i=(x+1):num_events
    beginning = ind(i) - window_around_event/2;
    last = ind(i)+ window_around_event/2;
    if beginning <1 || last>AmountOfData
        continue; %since can't add different sized window, just ignore
    end
    new_count=new_count+1;
    timepts = beginning:last;
    lastSpec=DataAfterCAR(:,freqs2plot, timepts);
    sum_windows = sum_windows + DataAfterCAR(:,freqs2plot, timepts);
    num_samples = num_samples + 1;
    last_used_ind = ind(i);
    current_num_events=prev_num_events+new_count;
    tmp=(lastSpec-averages)./stdevs;
    allStacked(:,current_num_events,:)=(lastSpec(:,6,:)-averages(:,6,:))./stdevs(:,6,:);
    %allStacked(:,:,:,current_num_events)=tmp;
    new_plot_flag=1;
end;
if num_samples>0
    num_samples = num_samples+old_num_samples;
    average = (old_average*old_num_samples + sum_windows)/(num_samples);
else %if no new events
    average = old_average;
    num_samples = old_num_samples;
    new_plot_flag=0;
end
%}
%%
% --- Executes on button press in loadsimdata.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to loadsimdata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of loadsimdata
load('C:\Users\Angela_2\Desktop\ANGELA\TDT_DataMAT\to_plot.mat')
load('C:\Users\Angela_2\Desktop\ANGELA\TDT_DataMAT\EC9_B58_start_ind_AND_threshold.mat')

load('C:\Users\Angela_2\Desktop\ANGELA\TDT_DataMAT\DataMAT.mat')
load('C:\Users\Angela_2\Desktop\ANGELA\baseline_stats_EC9_b58.mat')

vars=whos;
for i=1:length(vars)
    handles.(vars(i).name)=eval(vars(i).name);
end
printf('Done\n')
guidata(hObject, handles);


%%
%Misc functions created by GUIDE
function get_total_electrodes_Callback(hObject, eventdata, handles)
% hObject    handle to get_total_electrodes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns get_total_electrodes contents as cell array
%        contents{get(hObject,'Value')} returns selected item from get_total_electrodes
handles.number_of_electrodes_total = cellstr(get(hObject,'String'))
guidata(hObject, handles);
function get_total_electrodes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to get_total_electrodes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function time2collectBaseline_Callback(hObject, eventdata, handles)
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
function desired_freq_plot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to desired_freq_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function calculated_baseline_flag_Callback(hObject, eventdata, handles)
function threshold_Callback(hObject, eventdata, handles)
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
function desired_ANchan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to desired_ANchan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function popupmenu1_Callback(hObject, eventdata, handles)
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function numchannels_Callback(hObject, eventdata, handles)
function numchannels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numchannels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function number_of_electrodes_total_Callback(hObject, eventdata, handles)
function number_of_electrodes_total_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function number_of_electrodes_total_ButtonDownFcn(hObject, eventdata, handles)
function update_KeyPressFcn(hObject, eventdata, handles)
function eventCounter_ButtonDownFcn(hObject, eventdata, handles)


% --- Executes on selection change in get_num_analog.
function get_num_analog_Callback(hObject, eventdata, handles)
% hObject    handle to get_num_analog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns get_num_analog contents as cell array
%        contents{get(hObject,'Value')} returns selected item from get_num_analog


% --- Executes during object creation, after setting all properties.
function get_num_analog_CreateFcn(hObject, eventdata, handles)
% hObject    handle to get_num_analog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function desired_ANchan2_Callback(hObject, eventdata, handles)
% hObject    handle to desired_ANchan2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of desired_ANf desired_ANchan2 as a double


% --- Executes during object creation, after setting all properties.
function desired_ANchan2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to desired_ANchan2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function threshold2_Callback(hObject, eventdata, handles)
% hObject    handle to threshold2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of threshold2 as text
%        str2double(get(hObject,'String')) returns contents of threshold2 as a double


% --- Executes during object creation, after setting all properties.
function threshold2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threshold2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
