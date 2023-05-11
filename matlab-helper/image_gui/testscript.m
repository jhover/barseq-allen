cd 'C:\Program Files\Micro-Manager-2.0gamma'
%
clearvars
clc


import org.micromanager.*;
import mmcorej.*

mmc=CMMCore;
mmc.loadSystemConfiguration ('C:\Program Files\Micro-Manager-2.0gamma\barseqconfocal_Xc2-2.cfg');

cd D:\micromanager-matlab

%%
gui = StartMMStudio() 








%%
gui.show;
mmc = gui.getCore;
acq = gui.getAcquisitionEngine;
%%
clearvars
clc
import mmcorej.*;
mmc = CMMCore;
info = mmc.getVersionInfo();
info
%%
mmc.loadDevice("Camera","DemoCamera","DCam")

%%
mmc.initializeDevice("Camera")

%%
mmc.setCameraDevice("Camera");

%%

mmc.setExposure(50);
mmc.snapImage();
im=mmc.getImage();
%%
mmc.loadDevice("DHub","DemoCamera","DHub")
mmc.loadDevice("Camera","DemoCamera","DCam");
mmc.initializeDevice("DHub");

mmc.initializeDevice("Camera");


devices=mmc.getLoadedDevices();

%camera = mmc.getCameraDevice();
%props = mmc.getDevicePropertyNames(camera)

%%
mmc.loadSystemConfiguration ('C:\Program Files\Micro-Manager-2.0gamma\barseqconfocal_Xc2-2.cfg');