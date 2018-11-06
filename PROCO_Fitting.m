% PROCO MODEL DEFINITION
%**************************************************************************
% The Ohio State University
% Written by:   Michael Bush 
% Email:        bush.326@osu.edu, michael.bush1313@gmail.com
% Last update:  11/6/2018
%**************************************************************************
% This code reads a combined file of orthogonal training images and navigator positions, registers
% the training images, extracts translations and fits these translations
% with respect to navigator positions. It then creates a text file of model
% coefficients that can be read by a modified pulse sequence, allowing for
% real time updates of the imaging plane based on subject specific motion
% models. 

% The initial code segment requires proprietary software to perform image
% registration; we have provided a pre-run workspace for this segment to
% allow creation and visualization of the motion model for a single
% volunteer. 


% close all
% clear all
% clc
% addpath('moco_RA');
% 
% mkdir 'C:\Users\bush46\Documents\MATLAB\PROCOdata\MB2ch'
% mkdir 'C:\Users\bush46\Documents\MATLAB\PROCOdata\MBsax'
% 
% %This segment loads all image training data
% [im,orient,pixel_sizex,pixel_sizey,info,full_paths] = LoadDiccom('C:\Users\bush46\Documents\MATLAB\PROCOComboData\Combo');
% for i = 1:2:99
%     str_2ch = transpose(full_paths(:,i));
%     str_sax = transpose(full_paths(:,i+1));
%     copyfile(str_2ch,'C:\Users\bush46\Documents\MATLAB\PROCOdata\MB2ch');
%     copyfile(str_sax,'C:\Users\bush46\Documents\MATLAB\PROCOdata\MBsax');
% end
% 
% %Images are separated into 2ch and sax planes, pixel sizes and image
% %orientations are likewise extracted
% [im2ch,orient2ch,pixel_sizex_2ch,pixel_sizey_2ch] = LoadDiccom('C:\Users\bush46\Documents\MATLAB\PROCOdata\MB2ch');
% [imsax,orientsax,pixel_sizex_sax,pixel_sizey_sax] = LoadDiccom('C:\Users\bush46\Documents\MATLAB\PROCOdata\MBsax');
% u = [orient2ch(4),orient2ch(5),orient2ch(6)];
% v = [orientsax(4),orientsax(5),orientsax(6)];
% w = [orientsax(1),orientsax(2),orientsax(3)];
% SharedAxisAngle1 = atan2d(norm(cross(u,v)),dot(u,v))
% SharedAxisAngle2 = atan2d(norm(cross(u,w)),dot(u,w))
% B2ch = 10*im2ch/max(im2ch(:));   %histeq?
% Bsax = 10*imsax/max(imsax(:));
% 
% %Determination of orthogonal direction
% if abs(SharedAxisAngle2) > 1
%     orth_vector = 1
% else
%     orth_vector = 2
% end
% 
% %Reading in Navigator values
% fileID = fopen('C:\Users\bush46\Documents\MATLAB\PROCOComboData\navdata.txt');
% formatSpec = '%f %d';
% 
%     
% sizeA = [2 Inf];
% A = transpose(fscanf(fileID,formatSpec,sizeA));
% counter = 1;
% for k = 1:length(A(:,1))
%     if A(k,2) == length(im2ch(1,1,:))*2  
%         navdata(counter) = A(k,1);
%         counter = counter+1;
%     end
% end
% 
% %separating the continuous navigator signal into the correct temporal
% %location
% counter = 1;
% for i = 1:2:length(navdata)-1
%     nav2ch(counter) = navdata(i);
%     navsax(counter) = navdata(i+1);
%     counter = counter+1;
% end
% 
% %Selection of maximum expiratory position as the reference to register to
% [M2ch,ref2ch] = max(nav2ch);
% [Msax,refsax] = max(navsax);
% 
% close all
% 
% %Perform image registration here, extract deformation field parameters Dx
% %and Dy 
% param.mocoReg = 12;
% [Dx_2ch,Dy_2ch, DxInv_2ch, DyInv_2ch] = mocoDisp(B2ch, ref2ch, param.mocoReg);
% Am_2ch  = mocoApply(B2ch, Dx_2ch, Dy_2ch);
% [Dx_sax,Dy_sax, DxInv_sax, DyInv_sax] = mocoDisp(Bsax, refsax, param.mocoReg);
% Am_sax  = mocoApply(Bsax, Dx_sax, Dy_sax);
% x_xreg_2ch = cat(2, B2ch, Am_2ch);
% x_xreg_sax = cat(2, Bsax, Am_sax);

%Setup initial location for ROI
[Icrop_sax, rect_sax] = imcrop(imagesc(abs(squeeze(Bsax(:,:,refsax))),[0, max(abs(Bsax(:)))]));
[Icrop_2ch, rect_2ch] = imcrop(imagesc(abs(squeeze(B2ch(:,:,refsax))),[0, max(abs(B2ch(:)))]));

%sax analysis%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Loops through the deformation fields and extracts average translations
%within the ROI. 
for j = 1:length(Bsax(1,1,:))
    [Dx_crop_sax(:,:,j)] = pixel_sizex_sax*imcrop(Dx_sax(:,:,j), rect_sax);
    Dx_avg_sax(j) = mean(mean(Dx_crop_sax(:,:,j)));
    Dx_max_sax(j) = max(max(abs(Dx_crop_sax(:,:,j))));
    [Dy_crop_sax(:,:,j)] = pixel_sizey_sax*imcrop(Dy_sax(:,:,j), rect_sax);
    Dy_avg_sax(j) = mean(mean(Dy_crop_sax(:,:,j)));
    Dy_max_sax(j) = max(max(abs(Dy_crop_sax(:,:,j))));
end

%separation into respiratory phase based on direction of navigator curve
insp_count_sax = 1;
exp_count_sax = 1;
for k = 2:length(navsax)
    if navsax(k)-navsax(k-1)>0
        nav_exp_sax(exp_count_sax) = navsax(k); 
        exp_dx_sax(exp_count_sax) = Dx_avg_sax(k);
        exp_dy_sax(exp_count_sax) = Dy_avg_sax(k);
        exp_count_sax = exp_count_sax +1;
    else
        nav_insp_sax(insp_count_sax) = navsax(k); 
        insp_dx_sax(insp_count_sax) = Dx_avg_sax(k);
        insp_dy_sax(insp_count_sax) = Dy_avg_sax(k);
        insp_count_sax = insp_count_sax +1;
    end
end

%2ch analysis%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Repeat for 2ch frames
for j = 1:length(B2ch(1,1,:))
    [Dx_crop_2ch(:,:,j)] = pixel_sizex_2ch*imcrop(Dx_2ch(:,:,j), rect_2ch);
    Dx_avg_2ch(j) = mean(mean(Dx_crop_2ch(:,:,j)));
    Dx_max_2ch(j) = max(max(abs(Dx_crop_2ch(:,:,j))));
    [Dy_crop_2ch(:,:,j)] = pixel_sizey_2ch*imcrop(Dy_2ch(:,:,j), rect_2ch);
    Dy_avg_2ch(j) = mean(mean(Dy_crop_2ch(:,:,j)));
    Dy_max_2ch(j) = max(max(abs(Dy_crop_2ch(:,:,j))));
end

%separation into respiratory phase
insp_count_2ch = 1;
exp_count_2ch = 1;
for k = 2:length(nav2ch)
    if nav2ch(k)-nav2ch(k-1)>0
        nav_exp_2ch(exp_count_2ch) = nav2ch(k); 
        exp_dx_2ch(exp_count_2ch) = Dx_avg_2ch(k);
        exp_dy_2ch(exp_count_2ch) = Dy_avg_2ch(k);
        exp_count_2ch = exp_count_2ch +1;
    else
        nav_insp_2ch(insp_count_2ch) = nav2ch(k); 
        insp_dx_2ch(insp_count_2ch) = Dx_avg_2ch(k);
        insp_dy_2ch(insp_count_2ch) = Dy_avg_2ch(k);
        insp_count_2ch = insp_count_2ch +1;
    end
end

%%%%combining data into x y z image domain
%All Xi information comes from sax view
X_image_exp = exp_dy_sax;
X_exp_nav = nav_exp_sax;
X_image_insp = insp_dy_sax;
X_insp_nav = nav_insp_sax;

if orth_vector == 1
    %All Zi information comes from 2ch view; if orth_vector = 1, dy from
    %2ch view is image Z, dx from 2ch view is image Y. 
    Z_image_exp = exp_dy_2ch;
    Z_exp_nav = nav_exp_2ch;
    Z_image_insp = insp_dy_2ch;
    Z_insp_nav = nav_insp_2ch;
    
    %Yi information comes from both sax and 2ch view, doubling our samples
    %for Yi dimension. This section combines nav positions and translations
    %into a single vector.
    Y_image_exp(1:length(nav_exp_sax)) = exp_dx_sax;
    Y_image_exp(length(nav_exp_sax)+1:length(nav_exp_sax)+length(nav_exp_2ch)) = exp_dx_2ch;
    Y_exp_nav(1:length(nav_exp_sax)) = nav_exp_sax;
    Y_exp_nav(length(nav_exp_sax)+1:length(nav_exp_sax)+length(nav_exp_2ch)) = nav_exp_2ch;
    Y_image_insp(1:length(nav_insp_sax)) = insp_dx_sax;
    Y_image_insp(length(nav_insp_sax)+1:length(nav_insp_sax)+length(nav_insp_2ch)) = insp_dx_2ch;
    Y_insp_nav(1:length(nav_insp_sax)) = nav_insp_sax;
    Y_insp_nav(length(nav_insp_sax)+1:length(nav_insp_sax)+length(nav_insp_2ch)) = nav_insp_2ch;
else
    %dx from 2ch view is image Z, dy from 2ch view is image Y. 
    Z_image_exp = exp_dx_2ch;
    Z_exp_nav = nav_exp_2ch;
    Z_image_insp = insp_dx_2ch;
    Z_insp_nav = nav_insp_2ch;
    
    Y_image_exp(1:length(nav_exp_sax)) = exp_dx_sax;
    Y_image_exp(length(nav_exp_sax)+1:length(nav_exp_sax)+length(nav_exp_2ch)) = exp_dy_2ch;
    Y_exp_nav(1:length(nav_exp_sax)) = nav_exp_sax;
    Y_exp_nav(length(nav_exp_sax)+1:length(nav_exp_sax)+length(nav_exp_2ch)) = nav_exp_2ch;
    Y_image_insp(1:length(nav_insp_sax)) = insp_dx_sax;
    Y_image_insp(length(nav_insp_sax)+1:length(nav_insp_sax)+length(nav_insp_2ch)) = insp_dy_2ch;
    Y_insp_nav(1:length(nav_insp_sax)) = nav_insp_sax;
    Y_insp_nav(length(nav_insp_sax)+1:length(nav_insp_sax)+length(nav_insp_2ch)) = nav_insp_2ch;
end

%Perform fractional polynomial fitting for all models

%Xi exp 
for j = 1:length(X_exp_nav)
            QX_exp_nav(j,1) = 1;
            QX_exp_nav(j,2) = sqrt(X_exp_nav(j));
            QX_exp_nav(j,3) = X_exp_nav(j);
            QX_exp_nav(j,4) = 0;
end
a_exp_X = QX_exp_nav\transpose(X_image_exp);

%Xi insp
for j = 1:length(X_insp_nav)
            QX_insp_nav(j,1) = 1;
            QX_insp_nav(j,2) = sqrt(X_insp_nav(j));
            QX_insp_nav(j,3) = X_insp_nav(j);
            QX_insp_nav(j,4) = 0;
end
a_insp_X = QX_insp_nav\transpose(X_image_insp);

%Yi exp
for j = 1:length(Y_exp_nav)
            QY_exp_nav(j,1) = 1;
            QY_exp_nav(j,2) = sqrt(Y_exp_nav(j));
            QY_exp_nav(j,3) = Y_exp_nav(j);
            QY_exp_nav(j,4) = 0;
end
a_exp_Y = QY_exp_nav\transpose(Y_image_exp);

%Yi insp
for j = 1:length(Y_insp_nav)
            QY_insp_nav(j,1) = 1;
            QY_insp_nav(j,2) = sqrt(Y_insp_nav(j));
            QY_insp_nav(j,3) = Y_insp_nav(j);
            QY_insp_nav(j,4) = 0;
end
a_insp_Y = QY_insp_nav\transpose(Y_image_insp);

%Zi exp
for j = 1:length(Z_exp_nav)
            QZ_exp_nav(j,1) = 1;
            QZ_exp_nav(j,2) = sqrt(Z_exp_nav(j));
            QZ_exp_nav(j,3) = Z_exp_nav(j);
            QZ_exp_nav(j,4) = 0;
end
a_exp_Z = QZ_exp_nav\transpose(Z_image_exp);

%Zi insp
for j = 1:length(Z_insp_nav)
            QZ_insp_nav(j,1) = 1;
            QZ_insp_nav(j,2) = sqrt(Z_insp_nav(j));
            QZ_insp_nav(j,3) = Z_insp_nav(j);
            QZ_insp_nav(j,4) = 0;
end
a_insp_Z = QZ_insp_nav\transpose(Z_image_insp);

%This segment creates data to be plotted for model visualization. Each
%model will only be plotted over the range of binned and trained nav values for
%inspiration and expiration. 
nav_insp_mat = min([nav_insp_sax nav_insp_2ch]):0.01:max([nav_insp_sax nav_insp_2ch]);
nav_exp_mat = min([nav_exp_sax nav_exp_2ch]):0.01:max([nav_exp_sax nav_exp_2ch]);
for k = 1:length(nav_insp_mat)
    x_image_insp_model(k) = a_insp_X(1)+a_insp_X(2)*sqrt(nav_insp_mat(k))+a_insp_X(3)*nav_insp_mat(k);
    y_image_insp_model(k) = a_insp_Y(1)+a_insp_Y(2)*sqrt(nav_insp_mat(k))+a_insp_Y(3)*nav_insp_mat(k);
    z_image_insp_model(k) = a_insp_Z(1)+a_insp_Z(2)*sqrt(nav_insp_mat(k))+a_insp_Z(3)*nav_insp_mat(k);
end
for k = 1:length(nav_exp_mat)
    x_image_exp_model(k) = a_exp_X(1)+a_exp_X(2)*sqrt(nav_exp_mat(k))+a_exp_X(3)*nav_exp_mat(k);
    y_image_exp_model(k) = a_exp_Y(1)+a_exp_Y(2)*sqrt(nav_exp_mat(k))+a_exp_Y(3)*nav_exp_mat(k);
    z_image_exp_model(k) = a_exp_Z(1)+a_exp_Z(2)*sqrt(nav_exp_mat(k))+a_exp_Z(3)*nav_exp_mat(k);
end

%This segment transforms the model coefficients into the MRI physical
%coordinate system based on the orientation vectors provided in the DICOM
%headers of the sax and 2ch images. 
if orth_vector == 1
    a_exp_xphys = orientsax(1).*a_exp_X + orientsax(4).*a_exp_Y + orient2ch(1).*a_exp_Z;
    a_insp_xphys = orientsax(1).*a_insp_X + orientsax(4).*a_insp_Y + orient2ch(1).*a_insp_Z;
    a_exp_yphys = orientsax(2).*a_exp_X + orientsax(5).*a_exp_Y + orient2ch(2).*a_exp_Z;
    a_insp_yphys = orientsax(2).*a_insp_X + orientsax(5).*a_insp_Y + orient2ch(2).*a_insp_Z;
    a_exp_zphys = orientsax(3).*a_exp_X + orientsax(6).*a_exp_Y + orient2ch(3).*a_exp_Z;
    a_insp_zphys = orientsax(3).*a_insp_X + orientsax(6).*a_insp_Y + orient2ch(3).*a_insp_Z;
else
    a_exp_xphys = orientsax(1).*a_exp_X + orientsax(4).*a_exp_Y + orient2ch(4).*a_exp_Z;
    a_insp_xphys = orientsax(1).*a_insp_X + orientsax(4).*a_insp_Y + orient2ch(4).*a_insp_Z;
    a_exp_yphys = orientsax(2).*a_exp_X + orientsax(5).*a_exp_Y + orient2ch(5).*a_exp_Z;
    a_insp_yphys = orientsax(2).*a_insp_X + orientsax(5).*a_insp_Y + orient2ch(5).*a_insp_Z;
    a_exp_zphys = orientsax(3).*a_exp_X + orientsax(6).*a_exp_Y + orient2ch(6).*a_exp_Z;
    a_insp_zphys = orientsax(3).*a_insp_X + orientsax(6).*a_insp_Y + orient2ch(6).*a_insp_Z; 
end

%Plotting of data in physical domain for model visualization
for k = 1:length(nav_insp_mat)
    x_phys_insp_model(k) = a_insp_xphys(1)+a_insp_xphys(2)*sqrt(nav_insp_mat(k))+a_insp_xphys(3)*nav_insp_mat(k);
    y_phys_insp_model(k) = a_insp_yphys(1)+a_insp_yphys(2)*sqrt(nav_insp_mat(k))+a_insp_yphys(3)*nav_insp_mat(k);
    z_phys_insp_model(k) = a_insp_zphys(1)+a_insp_zphys(2)*sqrt(nav_insp_mat(k))+a_insp_zphys(3)*nav_insp_mat(k);
end
for k = 1:length(nav_exp_mat)
    x_phys_exp_model(k) = a_exp_xphys(1)+a_exp_xphys(2)*sqrt(nav_exp_mat(k))+a_exp_xphys(3)*nav_exp_mat(k);
    y_phys_exp_model(k) = a_exp_yphys(1)+a_exp_yphys(2)*sqrt(nav_exp_mat(k))+a_exp_yphys(3)*nav_exp_mat(k);
    z_phys_exp_model(k) = a_exp_zphys(1)+a_exp_zphys(2)*sqrt(nav_exp_mat(k))+a_exp_zphys(3)*nav_exp_mat(k);
end

figure(1)
plot(nav_insp_mat,x_image_insp_model,'r')
hold on
plot(nav_exp_mat,x_image_exp_model,'r--')
scatter(X_exp_nav,X_image_exp,'x','r')
scatter(X_insp_nav,X_image_insp,'o','r')

plot(nav_insp_mat,z_image_insp_model,'b')
plot(nav_exp_mat,z_image_exp_model,'b--')
scatter(Z_exp_nav,Z_image_exp,'x','b')
scatter(Z_insp_nav,Z_image_insp,'o','b')

plot(nav_insp_mat,y_image_insp_model,'k')
plot(nav_exp_mat,y_image_exp_model,'k--')
scatter(Y_exp_nav,Y_image_exp,'x','k')
scatter(Y_insp_nav,Y_image_insp,'o','k')
hold off
xlabel('Nav Position (mm)')
ylabel('Translation (mm)')
title('Fitting in image domain with combined Image Y data')

figure(2)
plot(nav_insp_mat,x_phys_insp_model,'r')
hold on
plot(nav_exp_mat,x_phys_exp_model,'r--')
plot(nav_insp_mat,y_phys_insp_model,'b')
plot(nav_exp_mat,y_phys_exp_model,'b--')
plot(nav_insp_mat,z_phys_insp_model,'k')
plot(nav_exp_mat,z_phys_exp_model,'k--')
hold off
title('Physical Domain Models')
legend('X Inspiration','X Expiration','Y Inspiration','Y Expiration','Z Inspiration','Z Expiration')
xlabel('Nav Position (mm)')
ylabel('Translation (mm)')

%Creation of the coefficient text file for transfer to the modified pulse
%sequence. Modifications must be made to the pulse sequence to read and
%store these coefficients for real-time imaging. 
c_imfit(2,1) = a_exp_xphys(1);
c_imfit(2,2) = a_exp_xphys(3);
c_imfit(2,3) = a_exp_xphys(4);
c_imfit(2,4) = a_exp_xphys(2);

c_imfit(2,5) = a_insp_xphys(1);
c_imfit(2,6) = a_insp_xphys(3);
c_imfit(2,7) = a_insp_xphys(4);
c_imfit(2,8) = a_insp_xphys(2);

c_imfit(2,9) = a_exp_yphys(1);
c_imfit(2,10) = a_exp_yphys(3);
c_imfit(2,11) = a_exp_yphys(4);
c_imfit(2,12) = a_exp_yphys(2);

c_imfit(2,13) = a_insp_yphys(1);
c_imfit(2,14) = a_insp_yphys(3);
c_imfit(2,15) = a_insp_yphys(4);
c_imfit(2,16) = a_insp_yphys(2);

c_imfit(2,17) = a_exp_zphys(1);
c_imfit(2,18) = a_exp_zphys(3);
c_imfit(2,19) = a_exp_zphys(4);
c_imfit(2,20) = a_exp_zphys(2);

c_imfit(2,21) = a_insp_zphys(1);
c_imfit(2,22) = a_insp_zphys(3);
c_imfit(2,23) = a_insp_zphys(4);
c_imfit(2,24) = a_insp_zphys(2);

fileID = fopen('motion_correction_coeffs.txt','w');
fprintf(fileID,'%6.2f   %12.8f\r\n',c_imfit);
fclose(fileID);
type motion_correction_coeffs_imfit.txt
