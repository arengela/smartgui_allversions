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

% Last Modified by GUIDE v2.5 01-Nov-2011 17:44:42

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

%{
start_flag = get(hObject,'value');
save start.mat
pause(.01)
%}

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

% Flags
event_flag = 1;  % Average plot
single_event_flag = 1;
single_stacked_flag = 1;
calculated_baseline_flag = get(handles.calculated_baseline_flag,'Value');

% Retrieve from GUI
time2collectBaseline = str2num(get(handles.time2collectBaseline,'String'));
desired_freq_plot = str2num(get(handles.desired_freq_plot,'String'));
freq_band_singlestacked = str2num(get(handles.freq_band_singlestacked,'String'));
desired_ANchan = str2num(get(handles.desired_ANchan,'String'));
threshold = str2num(get(handles.threshold,'String'));
subj_id = get(handles.subj_id,'String');
elec_tmp = get(handles.get_total_electrodes,'String');
number_of_electrodes_total=str2num(elec_tmp{get(handles.get_total_electrodes,'Value')});

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



pth = pwd;
addpath([pth filesep 'movie']);
rehash path;

%%
if calculated_baseline_flag==1
    %load baseline_stats
    %load('C:\Users\Connie\Dropbox\ChangLab\Users\Angela\SMART\baseline_stats_EC555_B2.mat');

    uiopen('load')
    medians = squeeze(medians(:,:,1));
end

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

DA = actxcontrol('TDevAcc.X');% Loads ActiveX methods
if DA.ConnectServer('Local')>0 %Checks to see if there is an OpenProject and Workbench load
    if DA.SetSysMode(2)==1%Starts OpenEx project running and acquiring data. 1 = Standby, 2 = Preview, !!MUST CHANGE TO 3 WHEN RECORDING
        if DA.GetSysMode>1 % Checks to see if the Project is running
            if not(isnan(DA.ReadTargetVEX(dproj_name{1},0,10,'F32','F64'))) %checks to see if it loaded the correct project
                for i=1:length(dproj_name)
                    BufferSize(i)=DA.GetTargetSize(dproj_name{i});% returns the Buffer size (bookkeeping)
                end
                if event_flag ==1 || single_event_flag==1
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
                    if event_flag ==1 || single_event_flag==1
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
                    NewDataMAT=zeros(number_of_channels*length(AcqPos),num_freq_bands,number_of_points2read/(num_freq_bands*number_of_channels)); ;
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
                    if event_flag ==1 || single_event_flag==1
                       % reshape ANdata by chan x time
                        ANNewDataMAT=reshape(ANNewData,ANnumber_of_channels, ANnumber_of_points2read/ANnumber_of_channels);
                        
                        % Place data into MATLAB buffers
                        % finalPos = posInNewData_AfterCAR+sampling_rate*number_of_sec2read-1;
                        finalPos = posInNewData_AfterCAR + number_of_points2read/(num_freq_bands*number_of_channels) - 1;
                        ANfinalPos = ANposInNewData + ANnumber_of_points2read/ANnumber_of_channels - 1;
                        
                        DataAfterCAR(:,:,posInNewData_AfterCAR:finalPos)=LogNewData;
                        ANNewData_finalMAT(:,ANposInNewData:ANfinalPos)=ANNewDataMAT;
                        
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
                            ind = find(trigger,detected_num_events);
                            ind(find([0 diff(ind)<(sampling_rate*1)]))=[]; %no events within 1 sec...
                            %of each other will be acknowledged
                            num_events = length(ind);
                            
                            set(0,'currentfigure',handles.figure1);
                            set(handles.figure1,'CurrentAxes',fh5)
                            hold on
                            plot(matlab_loopCounter,num_events,'r*');
                            
                            [eventRelatedAvg,num_samples,last_used_ind] = average_event_window(num_events, ind, window_around_event,...
                                DataAfterCAR,finalPos(1),number_of_electrodes_total,...
                                num_avgEvent_freqs2plot, avgEvent_freqs2plot,...
                                eventRelatedAvg, num_samples, last_used_ind);
                            
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
                        if detected_num_events>0 && single_event_flag == 1
                            ind = find(trigger,detected_num_events);
                            ind(find([0 diff(ind)<(sampling_rate*1)]))=[]; %no events within 1 sec...
                            %of each other will be acknowledged
                            num_events = length(ind);
                            
                            lastEvent = single_event_window(ind, window_around_event,...
                                DataAfterCAR,finalPos(1),number_of_electrodes_total,...
                                num_singleEvent_freqs2plot, singleEvent_freqs2plot);
                            ZscoreLastEvent=(lastEvent-averages)./stdevs;
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
                        end
%%
%plot single stacked
                        
                        if detected_num_events>0 && single_stacked_flag ==1
                            ind = find(trigger,detected_num_events);
                            ind(find([0 diff(ind)<(sampling_rate*1)]))=[]; %no events within 1 sec...
                            %of each other will be acknowledged
                            num_events = length(ind);
                            to_plot= SingleStackedTrials(DataAfterCAR, ind, num_events,...
                                averages, stdevs, freq_band_singlestacked, window_around_event);
                            num_events = size(to_plot,2);
                            if num_events ==0
                                continue
                            end
                            set(fh6,'DefaultAxesytick',[0.5: num_events:to_plot_grid*num_events])
                            set(0,'currentfigure',fh6);
                            flattened = reshape_3Ddata(to_plot, window_around_event,num_events, to_plot_grid);
                            to_plot = artifactrejection(flattened, window_around_event,-5.5, 5.5, .5, .8, to_plot_grid);
                            imagesc(real(to_plot),[-7 7]);
                            if num_events >1
                                fig_num_y_stacked = [1.1:num_events:to_plot_grid*num_events];
                                fig_num_y_stacked = reshape(repmat(fig_num_y_stacked, [to_plot_grid 1]), number_of_electrodes_total, 1);
                                text(fig_num_x, fig_num_y_stacked, fig_nums)
                                stackedfig_dash_y = repmat([0 num_events*to_plot_grid],to_plot_grid,1)';
                                line(fig_dash_x,stackedfig_dash_y,'LineStyle','--','linewidth',1','color',[0 0 0]);
                            end
                            grid on
                            
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
    end
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


% --- Executes on button press in Simulator.
function Simulator_Callback(hObject, eventdata, handles)
% hObject    handle to Simulator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Simulator

simulator_motormapping_SMART(handles)


% --- Executes on button press in loadsimdata.
function loadsimdata_Callback(hObject, eventdata, handles)
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
