function [data,nifti] = Bruker2nifti_smsEPI(pathData,scale,crop)
% Read Bruker scan info and fid to data struct
% and create nifti struct for save_untouch_nii()
%
% Input 
%       pathData  -  path to scan folder
%       scale     -  voxel size scaling factor, default = 10
%       crop      -  set to 1 to crop the image by MBfactor/2
%                           2 to attemp automatically detect volume        
%
% H.-L. Lee, 20190320
%

if nargin < 3
    crop = 0;
end

if nargin < 2
    scale = 10;
end

% Read data from scan folder
data = read_Bruker_raw(pathData);

if ~isfield(data,'kspace')
    data.kspace = ift2d(double(data.img));
    data.kspace = permute(data.kspace,[1 2 3 5 4]);
end

% Generate image array
if isfield(data.info,'MBfactor')    % Reconstruct SMS data
    data.kspace = sms_recon(data);
    imgo  = ft2d(data.kspace);
    if data.info.RefFlag == 1
        img            = abs(imgo(:,:,:,:,1));
        img(:,:,:,:,2) = abs(imgo(:,:,:,:,end));
    else
        img = imgo;
    end
    img = permute(sos(img,4),[1 2 3 5 4]);
    img = img*data.info.RECO_map_slope;
    % Crop image
    if crop == 1 && data.info.MBfactor>2
        Maty  = data.info.dims(2);
        ysize = Maty/data.info.MBfactor;
        img   = img(:,Maty/2-ysize+1:Maty/2+ysize,:,:,:);
        data.info.dims(2) = ysize*2;
        data.info.voxelOffset(2) = data.info.voxelOffset(2) + data.info.pixdim(2)*(Maty/2-ysize);
    elseif crop == 2 && data.info.MBfactor>2
        [y1,y2] = autocrop_y(img,data.info.dims(2)/data.info.MBfactor);
        img     = img(:,y1:y2,:,:,:); 
        data.info.dims(2) = y2 - y1 + 1;
        data.info.voxelOffset(2) = data.info.voxelOffset(2) + data.info.pixdim(2)*(y1-1);
    end
else                                % use 2dseq
    img = data.img;
end

% Change voxel size
data.info.idist       = data.info.idist*scale;
data.info.pixdim      = data.info.pixdim*scale;
data.info.offset      = data.info.offset*scale;
data.info.FOV         = data.info.FOV*scale;
data.info.voxelOffset = data.info.voxelOffset*scale;

% Fill in nifti header info
nifti.hdr.hist = fill_hdr_hist(data.info);
nifti.hdr.hk   = fill_hdr_hk(data.info);
nifti.hdr.dime = fill_hdr_dime(data.info);
nifti.filetype = 2;
nifti.untouch  = 1;

nifti.img                = int16(img);
nifti.hdr.dime.datatype  = 4;
nifti.hdr.dime.bitpix    = 16;
nifti.hdr.dime.dim(1)    = ndims(img);
nifti.hdr.dime.dim(2:3)  = size(img(:,:,1,1));
nifti.hdr.dime.dim(4)    = size(img,3);
nifti.hdr.dime.dim(5)    = size(img,4);
nifti.hdr.dime.glmax     = max(img(:));
nifti.hdr.dime.glmin     = min(img(:));
nifti.hdr.dime.scl_slope = 1/data.info.RECO_map_slope;

end


%% Local functions

%% Phase adjustment for SMS k-space
function ktemp = sms_recon(data)

rpEPI  = 1;    % 1 - fwd or rev EPI;  2 - interleaving fwd and rev EPI
yshift = 0.0;

Matx     = data.info.dims(1);
Maty     = data.info.dims(2);
Matz     = data.info.dims(3);
EncMaty  = data.info.PVM_EncMatrix(2);
rpEPI    = rpEPI * data.info.ACQ_scaling_phase;
Gzband   = data.info.MBfactor;
Gz       = data.info.RevGz;
refMode  = data.info.RefFlag;

if data.info.RFband == 0
    RFband = Gzband;
else
    RFband = 1;
end

FOVz           = data.info.idist*double(RFband)*double(Matz);
EpiReadCenter  = data.info.PVM_EpiReadCenter;

kspace = data.kspace;

kspace_regrid = kspace;
kspace_regrid = circshift(kspace_regrid,(Matx/2)-EpiReadCenter+1,1);

%  Correct for the extra dephasing caused by Gz blip gradient according to slice offset

halfFOVz     = 0.5*FOVz; 
sliceOffset = double(data.info.sliceOffsets);

ktemp  = zeros(Matx,Maty,Matz*RFband,size(kspace,4),size(kspace,5));
if refMode == 1     % Reference mode data
    if Gzband > 1
        for band = 1:RFband
            if RFband ~= 1
                sliceOffsetEff = sliceOffset - data.info.idist*length(sliceOffset)*(1-double(band)+(double(RFband)-1)/2);
            else
                sliceOffsetEff = sliceOffset;
            end
            for slice = 1:Matz
                for y = 1:EncMaty
                    ymod = double(mod(y-1,Gzband))-double((Gzband-1))/2;
                    ktemp(:,y+Maty-EncMaty,Matz*(band-1)+slice,:,1:2:end) = kspace_regrid(:,y+Maty-EncMaty,slice,:,1:2:end)*exp(1i*ymod*pi*sliceOffsetEff(slice)/halfFOVz);
                    ktemp(:,y+Maty-EncMaty,Matz*(band-1)+slice,:,2:2:end) = kspace_regrid(:,y+Maty-EncMaty,slice,:,2:2:end)*exp(-1i*ymod*pi*sliceOffsetEff(slice)/halfFOVz);
                end
            end
        end
    else
        ktemp = kspace_regrid;
    end
    ktemp = fovshift(ktemp,2,data.info.offset(2),data.info.FOV(2),1);
    ktemp(:,:,:,:,2:2:end) = rpEPI_fovshift(ktemp(:,:,:,:,2:2:end),2,data.info.offset(2), data.info.FOV(2), 1);
    ktemp = fovshift(ktemp,2,-yshift,Maty,1);
    ktemp = ktemp/4;
    
else                % Not reference mode
    if Gzband > 1
        for band = 1:RFband
            if RFband ~= 1
                sliceOffsetEff = sliceOffset - data.info.idist*length(sliceOffset)*(1-double(band)+(double(RFband)-1)/2);
            else
                sliceOffsetEff = sliceOffset;
            end
            for slice = 1:Matz
                for y = 1:EncMaty
                    ymod = double(mod(y-1,Gzband))-double(Gzband-1)/2;
                    ktemp(:,y+Maty-EncMaty,Matz*(band-1)+slice,:,:) = kspace_regrid(:,y+Maty-EncMaty,slice,:,:)*exp(Gz*1i*ymod*pi*sliceOffsetEff(slice)/halfFOVz);
                end
            end
        end
    else
        ktemp = kspace_regrid;
    end

    if rpEPI == -1
        ktemp = fovshift(ktemp,2,data.info.offset(2),data.info.FOV(2),1);
        ktemp = rpEPI_fovshift(ktemp,2,data.info.offset(2), data.info.FOV(2), 1);
        ktemp = fovshift(ktemp,2,-yshift,Maty,1);
    elseif rpEPI == 1
        ktemp = fovshift(ktemp,2,data.info.offset(2),data.info.FOV(2),1);
        ktemp = fovshift(ktemp,2,-yshift,Maty,1);
    elseif rpEPI == 2
        ktemp = fovshift(ktemp,2,data.info.offset(2),data.info.FOV(2),1);
        ktemp(:,:,:,:,2:2:end) = rpEPI_fovshift(ktemp(:,:,:,:,2:2:end),2,data.info.offset(2), data.info.FOV(2), 1);
        ktemp = fovshift(ktemp,2,-yshift,Maty,1);
    elseif rpEPI == -2
        ktemp = fovshift(ktemp,2,data.info.offset(2),data.info.FOV(2),1);
        ktemp(:,:,:,:,1:2:end) = rpEPI_fovshift(ktemp(:,:,:,:,1:2:end),2,data.info.offset(2), data.info.FOV(2), 1);
        ktemp = fovshift(ktemp,2,-yshift,Maty,1);
    end
end

end


%% Build nifti.hdr.dime struct
function dime = fill_hdr_dime(info)
    if info.dims(4) > 1
        dime.dim(1) = 4;
    else
        dime.dim(1) = 3;
    end
    dime.dim(2:5)    = info.dims;
    dime.dim(6:8)    = 1;
    dime.intent_p1   = 0;
    dime.intent_p2   = 0;
    dime.intent_p3   = 0;
    dime.intent_code = 0;
    dime.slice_start = 0;  
    dime.datatype    = info.datatype;
    dime.bitpix      = info.bitpix;
    dime.pixdim(1)   = 1;
    dime.pixdim(2:4) = info.pixdim;
    dime.pixdim(5)   = info.TR;
    dime.pixdim(6:8) = 0;
    dime.vox_offset  = 352;
    dime.scl_slope   = 1;
    dime.scl_inter   = 0;
    dime.slice_end   = 0;
    dime.slice_code  = info.PVM_ObjOrderScheme;
    dime.xyzt_units  = 10;
    dime.cal_max     = 0;
    dime.cal_min     = 0;
    dime.slice_duration = 0;
    dime.toffset        = 0;
end


%% Build nifti.hdr.hk struct
function hk = fill_hdr_hk(info)
    hk.sizeof_hdr    = 348;
    hk.data_type     = '';
    hk.db_name       = '';
    hk.extents       = 0;
    hk.session_error = 0;
    hk.regular       = 'r';
    hk.dim_info      = 0;
end


%% Build nifti.hdr.hist struct
function hist = fill_hdr_hist(info)
    hist.descrip     = info.descrip;
    hist.aux_file    = '';
    hist.qform_code  = 0;
    hist.sform_code  = 1;
    hist.quatern_b   = 0;
    hist.quatern_c   = 0;
    hist.quatern_d   = 0;
    hist.qoffset_x   = 0;
    hist.qoffset_y   = 0;
    hist.qoffset_z   = 0;
    hist.intent_name = '';
    hist.magic       = 'n+1';      
    
    temp = [info.GradOrient.*info.pixdim; info.voxelOffset'];
    hist.srow_x      =  temp(:,1)';
    hist.srow_y      =  temp(:,3)';
    hist.srow_z      = -temp(:,2)';      
end



%% Read scan info and fid
function data = read_Bruker_raw(pathName)
    % read scan info
    info = read_Bruker_info(pathName);
    
%     % read the original fid data
%     fid      = fopen([pathName '/fid'],'r');
%     fidTemp  = fread(fid,'int32','l');
%     data.fid = fidTemp(1:2:end) - 1i*fidTemp(2:2:end);
%     fclose(fid);

    % read 2dseq
    fid      = fopen([pathName '/pdata/1/2dseq'],'r');
    img      = fread(fid,'int16','l');
    data.img = reshape(img,info.dims(1:4));
    fclose(fid);

    % read fidREORD if exists
    if exist([pathName '/fidREORD0'],'file')
        data.kspace = zeros([info.dims(1:3) info.NReceivers info.dims(4)]);
        for c = 1:info.NReceivers
            fidName = [pathName '/fidREORD' num2str(c-1)];
            fid     = fopen(fidName,'r');
            fidTemp = fread(fid,'float64','l');
            fidCopy = fidTemp(1:2:end) - 1i*fidTemp(2:2:end);
            data.kspace(:,:,:,c,:) = reshape(fidCopy,[info.dims(1:3) 1 info.dims(4)]);
            fclose(fid);
        end
    end
    
    data.info = info;
end


%% Read scan info from acqp, method and reco files
function info = read_Bruker_info(pathName)
    
    [pathStudy,Enum,fileName] = fileparts(pathName);
    if numel(pathStudy) == 0
        pathE = Enum;
    else
        pathE = [pathStudy filesep Enum];
    end

    % Read acqp file
    pathACQP = [pathE filesep 'acqp'];
    fid      = fopen(pathACQP,'rt');

    while ~feof(fid)
        line    = fgetl(fid);
        parName = textscan(line,'%s','Delimiter','=');
        parName = parName{1};
        
        switch parName{1}
            case '##$ACQ_protocol_name'
                info.descrip = fgetl(fid);
                info.descrip = strip(info.descrip,'>');
                info.descrip = strip(info.descrip,'<');
            case '##$ACQ_dim'
                info.ACQ_dim = str2double(parName{2});
            case '##$ACQ_read_offset'
                line2          = strread(fgetl(fid),'%f'); 
                info.offset(1) = line2(1);           
            case '##$ACQ_phase1_offset' 
                line2          = strread(fgetl(fid),'%f'); 
                info.offset(2) = line2(1);  
            case '##$ACQ_slice_offset'
                sliceOffset = [];
                while true
                    line2 = fgetl(fid);
                    temp  = strread(line2,'%c','delimiter','#');
                    if temp(1)=='#' 
                        break; 
                    end
                    sliceOffset = [sliceOffset; strread(line2,'%f')];
                end
                info.sliceOffsets = sliceOffset;
                info.offset(3)    = mean(info.sliceOffsets);
            case '##$ACQ_slice_sepn'
                if info.ACQ_dim == 2
                    info.idist = eval(fgetl(fid));
                end
            case '##$NR' 
                info.dims(4) = eval(parName{2}); 
            case '##$ACQ_scaling_phase'
                info.ACQ_scaling_phase = eval(parName{2}); 
			case '##$NSLICES'
				info.dims(3) = eval(parName{2}); 
        end
    end
    fclose(fid);
    
    %Read method file
    pathMethod = [pathE filesep 'method'];
    fid        = fopen(pathMethod,'rt');
    
    while ~feof(fid)
        line    = fgetl(fid);
        parName = textscan(line,'%s','Delimiter','=');
        parName = parName{1};
        
        switch parName{1}
            case '##$PVM_RepetitionTime'        % TR in second
                info.TR = eval(parName{2})/1000;  
            case '##$PVM_ObjOrderScheme'        % slice order
                PVM_ObjOrderScheme = parName{2};  
                if strcmp(PVM_ObjOrderScheme,'Interlaced')
                    info.PVM_ObjOrderScheme = 3;
                elseif strcmp(PVM_ObjOrderScheme,'Sequential')
                    info.PVM_ObjOrderScheme = 1;
                end
            case '##$PVM_Fov'
                if strcmp(parName{2},'( 2 )')
                    [info.FOV(1), info.FOV(2)] = strread(fgetl(fid),'%f %f');                      
                elseif strcmp(parName{2},'( 3 )')
                    [info.FOV(1), info.FOV(2), info.FOV(3)] = strread(fgetl(fid),'%f %f %f');
                end
            case '##$PVM_SPackArrSliceOrient'
                info.sliceOrient  =   fgetl(fid);
            case '##$PVM_SPackArrReadOrient'
                info.readoutDir   =   fgetl(fid);
            case '##$PVM_SPackArrGradOrient'
                GradOrient = [];
                while true
                    line2  = fgetl(fid);
                    temp = strread(line2,'%c','delimiter','#');
                    if temp(1)=='#' 
                        break; 
                    end
                    GradOrient = [GradOrient;strread(deblank(line2))'];
                end
                info.GradOrient = reshape(GradOrient,3,3);
            case '##$PVM_EncNReceivers'
                info.NReceivers = eval(parName{2}); 
            case '##$PVM_EncMatrix'
                [mx,my] = strread(fgetl(fid),'%d %d');
                info.PVM_EncMatrix = [mx my];
            case '##$MBfactor'
                info.MBfactor = eval(parName{2}); 
            case '##$RevGz'
                info.RevGz  = eval(parName{2}); 
            case '##$RefFlag'
                RefFlag = parName{2};
                if strcmp(RefFlag,'Off')
                    info.RefFlag = 0;
                else
                    info.RefFlag = 1;
                end
            case '##$RFband'
                info.RFband = eval(parName{2});
            case '##$PVM_EncCentralStep1'
                info.PVM_EncCentralStep = eval(parName{2});
            case '##$PVM_EpiReadCenter'
                info.PVM_EpiReadCenter  = eval(parName{2});
        end
    end
    fclose(fid);
 
    %Read reco file
    pathReco = [pathE filesep 'pdata' filesep '1' filesep 'reco'];
    fid      = fopen(pathReco,'rt');
    
    while ~feof(fid) 
        line    = fgetl(fid);
        parName = textscan(line,'%s','Delimiter','=');
        parName = parName{1};
        
        switch parName{1}
            case '##$RECO_map_slope'            % image scaling factor
                line2 = fgetl(fid);
                slope = strread(line2,'%f');
                info.RECO_map_slope = slope(1); 
            case '##$RECO_size'                 % image dimensions
                if strcmp(parName{2},'( 2 )')
                    [info.dims(1), info.dims(2)] = strread(fgetl(fid),'%d %d');    
                elseif strcmp(parName{2},'( 3 )')
                    [info.dims(1), info.dims(2), info.dims(3)] = strread(fgetl(fid),'%d %d %d');
                end         
            case '##$RECO_wordtype'
                info.wordtype = parName{2};
                switch info.wordtype
                    case '_8BIT_USGN_INT'
                        info.wordtype = 'uint8';
                        info.datatype = 2;
                        info.bitpix   = 8;
                    case '_16BIT_SGN_INT' 
                        info.wordtype = 'int16'; 
                        info.datatype = 4;
                        info.bitpix   = 16;
                    case '_32BIT_SGN_INT' 
                        info.wordtype = 'int32'; 
                        info.datatype = 8;
                        info.bitpix   = 32;
                    case '_32BIT_FLT' 
                        info.wordtype = 'float32';
                        info.datatype = 16;
                        info.bitpix   = 32;
                    case '_64BIT_FLT' 
                        info.wordtype = 'float64'; 
                        info.datatype = 64;
                        info.bitpix   = 64;
                    case '_8BIT_SGN_INT' 
                        info.wordtype = 'int8'; 
                        info.datatype = 256;
                        info.bitpix   = 8;
                    case '_16BIT_USGN_INT' 
                        info.wordtype = 'uint16'; 
                        info.datatype = 512;
                        info.bitpix   = 16;
                    case '_32BIT_USGN_INT' 
                        info.wordtype = 'uint32';  
                        info.datatype = 768;
                        info.bitpix   = 32;            
                end
        end  
    end 
    fclose(fid);
    
    % Spatial resolution
    info.pixdim(1:2) = (info.FOV(1:2)./info.dims(1:2))';
    if info.ACQ_dim == 2
        info.FOV(3)    = info.dims(3)*info.idist;
        info.pixdim(3) = info.idist;
    elseif info.ACQ_dim == 3
        info.pixdim(3) = info.FOV(3)/info.dims(3);
        info.idist     = info.pixdim(3);
    end
    
    if isfield(info,'MBfactor') 
        info.FOV(3) = info.FOV(3)*info.MBfactor;
    end
    
    % Orientation matrix
    switch info.sliceOrient
        case 'axial' 
            if strcmp(info.readoutDir,'L_R')                     
                xyzOrder = [1,2,3];              
            elseif   strcmp(info.readoutDir,'A_P')
                xyzOrder = [2,1,3];               
            end
        case 'sagittal'
            if strcmp(info.readoutDir,'H_F')                 
                xyzOrder = [2,1,3];            
            elseif strcmp(info.readoutDir,'A_P')
                xyzOrder = [1,2,3];               
            end
        case 'coronal'
            if strcmp(info.readoutDir,'H_F')                 
                xyzOrder = [2,1,3];             
            elseif strcmp(info.readoutDir,'L_R')
                xyzOrder = [1,2,3];               
            end        
    end  
    GradOrientNew    = info.GradOrient(:,xyzOrder); 
    info.offset(3)   = -info.offset(3);
    voxelOffset      = -info.offset(xyzOrder) - info.FOV(xyzOrder)/2 + info.pixdim(xyzOrder)/2;
    info.voxelOffset = GradOrientNew*voxelOffset';
    info.GradOrient  = GradOrientNew;
    
end


function img = ft2d(kdata)
    kdata = ifftshift(ifftshift(kdata,1),2);
    img   = fftshift(fftshift(fft2(kdata),1),2)/size(kdata,1)/size(kdata,2);
end

function img = ift2d(kdata)
    kdata = ifftshift(ifftshift(kdata,1),2);
    img   = fftshift(fftshift(ifft2(kdata),1),2);
end

function img = rpEPI_fovshift(data, dim, d, fov, k_flag)
    if nargin < 5
        k_flag = 0;
    end

    mat_x    = size(data,dim);
    A_index  = -floor(mat_x/2):(ceil(mat_x/2) - 1);
    A_column = (exp(1i*double(A_index)*2*pi*double(d)*2/double(fov))).';

    dim_index    = 1:length(size(data));
    dim_index    = circshift(dim_index,(dim-1),2);
    rep_num      = size(data);
    rep_num(dim) = 1;
    rep_num      = circshift(rep_num,(1-dim),2);
    A_mat        = repmat(A_column,rep_num);
    A_mat        = permute(A_mat,dim_index);

    data = flip(data,dim);

    if k_flag == 0
        k   = fftshift(fft(data, [], dim),dim);
        k   = k.*A_mat;
        img = abs(ifft(ifftshift(k, dim),[], dim));
    else
        k = data;
        img   = k.*A_mat;
    end
end


function img = fovshift(data, dim, d, fov, k_flag)
    if nargin < 5
        k_flag = 0;
    end

    mat_x    = size(data,dim);
    A_index  = -floor(mat_x/2):(ceil(mat_x/2) - 1); 
    A_column = exp(1i*A_index*2*pi*d/double(fov));

    rep_num      = size(data);
    rep_num(dim) = 1;

    A_mat = repmat(A_column,rep_num);

    if k_flag == 0
        k = fftshift(fft(data, [], dim),dim);
    else
        k = data;
    end

    k = k.*A_mat;

    if k_flag == 0
        img = abs(ifft(ifftshift(k, dim),[], dim));
    else
        img = k;
    end
end


function sosimg = sos(x ,dim, pnorm)
    if nargin < 2
        dim = size(size(x),2);
    end

    if nargin < 3
        pnorm = 2;
    end

    sosimg = (sum(abs(x).^pnorm,dim)).^(1/pnorm);
end

function [y1,y2] = autocrop_y(img,ywin)
    img    = mean(img,4);
    [mx,my,mz] = size(img);
    T      = adaptthresh(int16(img),0.4);
    mask   = mask_kmeans(T);
    mask(T<median(T(T>median(T(T>median(T(:))))))) = 0;
    marker = zeros(size(img));
    marker(mx/2,my/2,mz/2) = 1;
    im     = imreconstruct(marker,mask);
    temp   = find(sum(sum(im,1),3)>0);
    y1     = temp(1) - 6; 
    y2     = temp(end) + 1;
    if y2 > my
        y2 = my;
    end
    if y2-y1 < ywin-1
        y1 = y2 - ywin + 1;
    end
    if y1<1
        y1 = 1;
    end
end

function mask = mask_kmeans(img)

    [idx,c] = kmeans(img(:),2,'MaxIter',100);

    mask = reshape(idx,size(img));


    if c(1)>c(2)
        mask(mask==2) = 0;
    else
        mask = floor(mask/2);
    end 

    se = strel('disk',3);
    mask = imerode(mask, se);
    se = strel('disk',2);
    mask = imdilate(mask, se);
end