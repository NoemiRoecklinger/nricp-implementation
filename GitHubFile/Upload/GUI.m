function varargout = GUI(varargin)
% GUI MATLAB code for GUI.fig
%      GUI, by itself, creates a new GUI or raises the existing
%      singleton*.
%
%      H = GUI returns the handle to a new GUI or the handle to
%      the existing singleton*.
%
%      GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI.M with the given input arguments.
%
%      GUI('Property','Value',...) creates a new GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI

% Last Modified by GUIDE v2.5 17-May-2018 20:25:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_OutputFcn, ...
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


% --- Executes just before GUI is made visible.
function GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI (see VARARGIN)

% Choose default command line output for GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Normals.
function Normals_Callback(hObject, eventdata, handles)
% hObject    handle to Normals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of Normals
    
    global useNormals;
    if(get(hObject, 'Value') == 0.0)
        
        useNormals = 0;
    else
        
        useNormals = 1;
    end
    assignin('base','useNormals', useNormals);

% --- Executes on selection change in Source.
function Source_Callback(hObject, eventdata, handles)
% hObject    handle to Source (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns Source contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Source

    global src
    contents = cellstr(get(hObject,'String'));
    src = contents{get(hObject, 'Value')};
    
    assignin('base','src', src);


% --- Executes during object creation, after setting all properties.
function Source_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Source (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Target.
function Target_Callback(hObject, eventdata, handles)
% hObject    handle to Target (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Target contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Target

    global trgt;
    contents = cellstr(get(hObject,'String'));
    trgt = contents{get(hObject, 'Value')};
    
    assignin('base','trgt', trgt);

% --- Executes during object creation, after setting all properties.
function Target_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Target (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in NRICP.
function NRICP_Callback(hObject, eventdata, handles)
% hObject    handle to NRICP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    class = MainClass;
    class.execute();
    handle = findobj('Tag', 'result');
    axes(handle)
    
    global residuals;
    plot([1:size(residuals,2)], residuals,'r'), 
    title('Residuals over all iterations'), 
    xlabel('Number of iterations'), ylabel('Residuals'),
    legend('||AX-B||F2')
    
    global normal_mean;
    normalLabel = sprintf('Normal Mean = %f', normal_mean);
    set(handles.normalMean, 'String', normalLabel);
    
    global curv_diff;
    curvatureLabel = sprintf('Curvature Difference = %f', curv_diff);
    set(handles.Curvature, 'String', curvatureLabel);

% --- Executes on button press in ignoreBoundary.
function ignoreBoundary_Callback(hObject, eventdata, handles)
% hObject    handle to ignoreBoundary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ignoreBoundary

    global ignoreBoundary;
    if(get(hObject, 'Value') == 0.0)
        
        ignoreBoundary = 0;
    else
        
        ignoreBoundary = 1;
    end
    assignin('base','ignoreBoundary', ignoreBoundary);



% --- Executes on button press in rigidInit.
function rigidInit_Callback(hObject, eventdata, handles)
% hObject    handle to rigidInit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rigidInit
    global rigidInit;
    if(get(hObject, 'Value') == 0.0)
        
        rigidInit = 0;
    else
        
        rigidInit = 1;
    end
    assignin('base','rigidInit', rigidInit);



% --- Executes on selection change in MaxAlpha.
function MaxAlpha_Callback(hObject, eventdata, handles)
% hObject    handle to MaxAlpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns MaxAlpha contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MaxAlpha

    global MaxAlpha;
    contents = cellstr(get(hObject,'String'));
    MaxAlpha = contents{get(hObject, 'Value')};
    MaxAlpha = str2double(MaxAlpha);
    assignin('base','MaxAlpha', MaxAlpha);


% --- Executes during object creation, after setting all properties.
function MaxAlpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxAlpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in MinAlpha.
function MinAlpha_Callback(hObject, eventdata, handles)
% hObject    handle to MinAlpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns MinAlpha contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MinAlpha

    global MinAlpha;
    contents = cellstr(get(hObject,'String'));
    MinAlpha = contents{get(hObject, 'Value')};
    MinAlpha = str2double(MinAlpha);
    
    assignin('base','MinAlpha', MinAlpha);


% --- Executes during object creation, after setting all properties.
function MinAlpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MinAlpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in StepsAlpha.
function StepsAlpha_Callback(hObject, eventdata, handles)
% hObject    handle to StepsAlpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns StepsAlpha contents as cell array
%        contents{get(hObject,'Value')} returns selected item from StepsAlpha

    global StepsAlpha;
    contents = cellstr(get(hObject,'String'));
    StepsAlpha = contents{get(hObject, 'Value')};
    StepsAlpha = str2double(StepsAlpha);
    
    assignin('base','StepsAlpha', StepsAlpha);

% --- Executes during object creation, after setting all properties.
function StepsAlpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StepsAlpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Epsilon_Callback(hObject, eventdata, handles)
% hObject    handle to Epsilon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global Epsilon;
    get(hObject,'String');
    Epsilon = str2double(get(hObject,'String'));
    
    assignin('base','Epsilon', Epsilon);


% --- Executes during object creation, after setting all properties.
function Epsilon_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Epsilon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    evalin('base', 'clear variables');
    close all
    clc
    run GUI


% --- Executes on button press in close.
function close_Callback(hObject, eventdata, handles)
% hObject    handle to close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    evalin('base', 'clear variables');
    close all
    clc
