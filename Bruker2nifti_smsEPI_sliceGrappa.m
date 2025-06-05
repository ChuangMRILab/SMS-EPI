function [data,nifti] = Bruker2nifti_smsEPI_sliceGrappa(pathData,pathRef,kSize,scale)
% Read Bruker scan info and fid to data struct
% and create nifti struct for save_untouch_nii()
%
% Input 
%       pathData  -  path to scan folder
%       pathRef   -  path to referenec scan folder
%       kSize     -  Parallel reconstruction kernel size (1x2 array)
%                    default: [11 11]
%       scale     -  voxel size scalar
%
% H.-L. Lee, 20190516
%
if nargin < 4
    scale = 10;
end

if nargin < 3
    kSize = [11 11];
end

% Read reference from scan folder
dataRef = read_bruker_raw(pathRef);
dataRef.kspace = dataRef.kspace(:,:,:,:,1);

% Read data from scan folder
data = read_bruker_raw(pathData);
if ~isfield(data,'kspace')
    data.kspace = ift2d(double(data.img));
end

% Generate image array
if isfield(data.info,'MBfactor')    % Reconstruct SMS data
    if data.info.dims(3)*data.info.MBfactor ~= dataRef.info.dims(3)
        disp('Error. Reference and data do not match.')
        return;
    end
    
    dataRef.info.MBfactor = data.info.MBfactor;
    kRef = generateRef(dataRef);
    if data.info.MBfactor > 1
        data.kspace = sms_recon(data);
    end
    data.kspace = sliceGRAPPA(data.kspace, kRef, kSize);
    
    imgo = ft2d(data.kspace);
    img  = squeeze(sos(imgo,4));
    img  = img*30000/max(img(:));
else                                % use 2dseq
    img = data.img;
end

% Change voxel size
data.info.idist  = data.info.idist*scale;
data.info.resol  = data.info.resol*scale;
data.info.offset = data.info.offset*scale;
data.info.FOV    = data.info.FOV*scale;
data.info.pos0   = data.info.pos0*scale;

% Fill in nifti header info
nifti.hdr.hist = fill_hdr_hist(data.info);
nifti.hdr.hk   = fill_hdr_hk(data.info);
nifti.hdr.dime = fill_hdr_dime(data.info);
nifti.filetype = 2;
nifti.untouch  = 1;

nifti.img                = int16(img);
nifti.hdr.dime.dim(1)    = ndims(img);
nifti.hdr.dime.dim(2:4)  = size(img(:,:,:,1));
nifti.hdr.dime.dim(5)    = size(img,4);
nifti.hdr.dime.glmax     = max(img(:));
nifti.hdr.dime.glmin     = min(img(:));
nifti.hdr.dime.scl_slope = data.info.RECO_map_maxima/max(img(:));

% Parse file name
[pathstr,Enum,fileName] = fileparts(pathData);
newFileName             = [data.info.descrip 'X' Enum '_sGreco.nii.gz'];

% Save nifti to file
save_untouch_nii(nifti,newFileName);

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
slice_offset = double(data.info.slice_offsets);

ktemp  = zeros(Matx,Maty,Matz*RFband,size(kspace,4),size(kspace,5));
if refMode == 1     % Reference mode data
    if Gzband > 1
        for band = 1:RFband
            if RFband ~= 1
                slice_offseteff = slice_offset - data.info.idist*length(slice_offset)*(1-double(band)+(double(RFband)-1)/2);
            else
                slice_offseteff = slice_offset;
            end
            for slice = 1:Matz
                for y = 1:EncMaty
                    ymod = double(mod(y-1,Gzband))-double((Gzband-1))/2;
                    ktemp(:,y+Maty-EncMaty,Matz*(band-1)+slice,:,1:2:end) = kspace_regrid(:,y+Maty-EncMaty,slice,:,1:2:end)*exp(1i*ymod*pi*slice_offseteff(slice)/halfFOVz);
                    ktemp(:,y+Maty-EncMaty,Matz*(band-1)+slice,:,2:2:end) = kspace_regrid(:,y+Maty-EncMaty,slice,:,2:2:end)*exp(-1i*ymod*pi*slice_offseteff(slice)/halfFOVz);
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
                slice_offseteff = slice_offset - data.info.idist*length(slice_offset)*(1-double(band)+(double(RFband)-1)/2);
            else
                slice_offseteff = slice_offset;
            end
            for slice = 1:Matz
                for y = 1:EncMaty
                    ymod = double(mod(y-1,Gzband))-double(Gzband-1)/2;
                    ktemp(:,y+Maty-EncMaty,Matz*(band-1)+slice,:,:) = kspace_regrid(:,y+Maty-EncMaty,slice,:,:)*exp(Gz*1i*ymod*pi*slice_offseteff(slice)/halfFOVz);
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


%% Create reference k-space from single-band data
function ktemp = generateRef(data)

rpEPI  = 1;    % 1 - fwd or rev EPI;  2 - interleaving fwd and rev EPI
yshift = 0.0;

Matx     = data.info.dims(1);
Maty     = data.info.dims(2);
Matz     = data.info.dims(3);
EncMaty  = data.info.PVM_EncMatrix(2);
rpEPI    = rpEPI * data.info.ACQ_scaling_phase;
Gzband   = data.info.MBfactor;
Gz       = data.info.RevGz;

FOVz           = data.info.idist*double(Matz);
EpiReadCenter  = data.info.PVM_EpiReadCenter;

kspace = data.kspace;

kspace_regrid = kspace;
kspace_regrid = circshift(kspace_regrid,(Matx/2)-EpiReadCenter+1,1);

%  Correct for the extra dephasing caused by Gz blip gradient according to slice offset

halfFOVz     = 0.5*FOVz; 
slice_offset = double(data.info.slice_offsets);
nSlice       = Matz/Gzband;

ktemp  = zeros(Matx,Maty,Matz,size(kspace,4),2);
% Not reference mode
for band = 1:Gzband
    slice = (band-1)*nSlice+1:band*nSlice;
    for n = 1:Gzband
        slice_offseteff((n-1)*nSlice+1:n*nSlice) = slice_offset(slice)-slice_offset((n-1)*nSlice+1:n*nSlice);
    end
    for n = 1:nSlice
        for y = 1:EncMaty
            ymod = mod(y-1,Gzband)-(Gzband-1)/2;
            ktemp(:,y+Maty-EncMaty,slice(n),:,1) = kspace_regrid(:,y+Maty-EncMaty,slice(n),:,1)*exp(Gz*1i*ymod*pi*slice_offseteff(slice(n))/halfFOVz);
            for k = 1:Gzband
                ktemp(:,y+Maty-EncMaty,slice(n),:,2) = ktemp(:,y+Maty-EncMaty,slice(n),:,2) + kspace_regrid(:,y+Maty-EncMaty,(k-1)*nSlice+n,:,1)*exp(Gz*1i*ymod*pi*slice_offseteff((k-1)*nSlice+n)/halfFOVz);
            end
            ktemp(:,y+Maty-EncMaty,slice(n),:,2) = ktemp(:,y+Maty-EncMaty,slice(n),:,2) - ktemp(:,y+Maty-EncMaty,slice(n),:,1);
        end
    end
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




%% Build nifti.hdr.dime struct
function dime = fill_hdr_dime(par)
    if par.dims(4) > 1
        dime.dim(1) = 4;
    else
        dime.dim(1) = 3;
    end
    dime.dim(2:5)    = par.dims;
    dime.dim(6:8)    = 1;
    dime.intent_p1   = 0;
    dime.intent_p2   = 0;
    dime.intent_p3   = 0;
    dime.intent_code = 0;
    dime.slice_start = 0;  
    dime.datatype    = par.datatype;
    dime.bitpix      = par.bitpix;
    dime.pixdim(1)   = 1;
    dime.pixdim(2:4) = par.resol;
    dime.pixdim(5)   = par.TR;
    dime.pixdim(6:8) = 0;
    dime.vox_offset  = 352;
    dime.scl_slope   = 1;
    dime.scl_inter   = 0;
    dime.slice_end   = 0;
    dime.slice_code  = par.PVM_ObjOrderScheme;
    dime.xyzt_units  = 10;
    dime.cal_max     = 0;
    dime.cal_min     = 0;
    dime.slice_duration = 0;
    dime.toffset        = 0;
end


%% Build nifti.hdr.hk struct
function hk = fill_hdr_hk(par)
    hk.sizeof_hdr    = 348;
    hk.data_type     = '';
    hk.db_name       = '';
    hk.extents       = 0;
    hk.session_error = 0;
    hk.regular       = 'r';
    hk.dim_info      = 0;
end


%% Build nifti.hdr.hist struct
function hist = fill_hdr_hist(par)
    hist.descrip     = par.descrip;
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
    
    temp = [par.m_or.*par.resol; par.pos0'];
    hist.srow_x      =  temp(:,1)';
    hist.srow_y      =  temp(:,3)';
    hist.srow_z      = -temp(:,2)';      
end



%% Read scan info and fid
function A = read_bruker_raw(PathName)
    pars = get_Bruker_pars(PathName);
    
    % read the original fid data
    if exist([PathName '/fid'])
        fpre = fopen([PathName '/fid'],'r');
        A.fid = fread(fpre,'int32','l');
        A.fid = A.fid(1:2:end) - 1i * A.fid(2:2:end);
        fclose(fpre);
    end

    % read 2dseq
    if exist([PathName '/pdata/1/2dseq'])
        fpre    = fopen([PathName '/pdata/1/2dseq'],'r');
        img     = fread(fpre,'int16','l');
        A.img   = reshape(img,pars.dims(1),pars.dims(2),pars.dims(3),pars.dims(4));
        fclose(fpre);
    end

    % read fidCopy if exist
    if exist([PathName '/fidREORD0'])
        A.kspace = zeros(pars.dims(1),pars.dims(2),pars.dims(3),pars.NReceivers,pars.dims(4));
        for c = 1:pars.NReceivers
            fidName = [PathName '/fidREORD' num2str(c-1)];
            fpre = fopen(fidName,'r');
            fid = fread(fpre,'float64','l');
            fidCopy = fid(1:2:end) - 1i * fid(2:2:end);
            k = reshape(fidCopy,pars.dims(1),pars.dims(2),pars.dims(3),1,pars.dims(4));
            A.kspace(:,:,:,c,:) = k;
            fclose(fpre);
        end
    end
    
    A.info = pars;
end


%% Read scan info from acqp, method and reco files
function pars = get_Bruker_pars(path)

    pars.acq_dim     =   0;
    pars.scale       =   ones(4,1);
    pars.cmpx        =   0;
    pars.n_coils     =   1;
    pars.orient      =   '';
    pars.r_out       =   '';
    pars.idist       =   0;
    pars.m_or        =   zeros(3,2);
    pars.dims        =   zeros(4,1);
    pars.FOV         =   zeros(3,1);   
    pars.resol       =   zeros(3,1);   
    pars.offset      =   zeros(3,1);  
    pars.tp          =   '';
    pars.endian      =   'l';
    pars.n_acq       =   '';
    pars.TR          =   0;
    pars.PVM_ObjOrderScheme = 0;
    
    [pathstr,Enum,fileName] = fileparts(path);
    if numel(pathstr) == 0
        pathscan    =   Enum;
    else
        pathscan    =   [pathstr filesep Enum];
    end
    
    path_acqp   =   deblank([pathscan filesep 'acqp']);
    fid         =   fopen(deblank(path_acqp), 'rt');
 
    %Read acqp
    while feof(fid) == 0
        line    =   fgetl(fid);
        tag     =   strread(line,'%s','delimiter','=');
        switch tag{1}
            case '##$ACQ_protocol_name'
                pars.descrip     = fgetl(fid);
                pars.descrip     = strip(pars.descrip,'>');
                pars.descrip     = strip(pars.descrip,'<');
            case '##$ACQ_dim'
                pars.acq_dim     =   str2num(tag{2});
            case '##$ACQ_slice_sepn'
                if pars.acq_dim==2
                    pars.idist       =   eval(fgetl(fid));
                    pars.resol(3)    =   pars.idist;
                end
%             case '##$ACQ_slice_thick' 
%                 if pars.acq_dim==2              
%                     pars.resol(3)    =   eval(tag{2});
%                 end
            case '##$ACQ_slice_offset'
                          slices        =   [];
                          while true
                            read        =   fgetl(fid);
                            header      =   strread(read,'%c','delimiter','#');
                            if strcmp(header(1),'#') ||strcmp(header(1),'$') break; end
                            slices   =   [slices;strread(read,'%f','delimiter','\\ ')];
                          end
                          pars.slice_offsets =   slices;
                          pars.offset(3)     =   mean(pars.slice_offsets);
            case '##$ACQ_read_offset'
               read              =   strread(fgetl(fid),'%f','delimiter','\\ '); 
               pars.offset(1)    =   read(1);           
            case '##$ACQ_phase1_offset' 
               read              =   strread(fgetl(fid),'%f','delimiter','\\ '); 
               pars.offset(2)    =   read(1);  
            case '##$NR' 
               pars.dims(4)      =   eval(tag{2}); 
            case '##$ACQ_scaling_phase'
                pars.ACQ_scaling_phase = eval(tag{2}); 
			case '##$NSLICES'
				pars.dims(3) =  eval(tag{2}); 
        end
    end
    fclose(fid);
    
    
    %Read method
    path_method    =   deblank([pathscan filesep 'method']);
    fid            =   fopen(path_method, 'rt');
    while feof(fid) == 0
        line    =   fgetl(fid);
        tag     =   strread(line,'%s','delimiter','=');
        switch tag{1}
            case '##$PVM_RepetitionTime'
                pars.TR      =   eval(tag{2})/1000; 
%             case '##$PVM_NRepetitions' 
%                 pars.dims(4)     =   eval(tag{2}); 
            case '##$PVM_ScanTime'
                scan_time    =   eval(tag{2}); 
            case '##$PVM_SPackArrSliceOrient'
                pars.orient  =   fgetl(fid);
            case '##$PVM_SPackArrReadOrient'
                pars.r_out   =   fgetl(fid);
            case '##$PVM_SPackArrGradOrient'
                m = 0;
                while true
                    read    =   fgetl(fid);
                    header  =   strread(read,'%c','delimiter','#');
                    if header(1)=='#' break; end
                    m       =   [m;strread(deblank(read))'];
                end
                pars.m_or   =   m(2:10);
                pars.m_or   =   reshape(pars.m_or,[3,3]);
            case '##$PVM_EncNReceivers'
                pars.NReceivers = eval(tag{2}); 
            case '##$PVM_EncMatrix'
                [a,b] = strread(fgetl(fid),'%d %d');
                pars.PVM_EncMatrix = [a b];
            case '##$MBfactor'
                pars.MBfactor = eval(tag{2}); 
            case '##$RevGz'
                pars.RevGz  = eval(tag{2}); 
            case '##$RefFlag'
                RefFlag = tag{2};
                if strcmp(RefFlag,'Off')
                    pars.RefFlag = 0;
                else
                    pars.RefFlag = 1;
                end
            case '##$RFband'
                pars.RFband = eval(tag{2});
            case '##$PVM_EncCentralStep1'
                pars.PVM_EncCentralStep = eval(tag{2});
            case '##$PVM_EpiReadCenter'
                pars.PVM_EpiReadCenter  = eval(tag{2});
            case '##$PVM_ObjOrderScheme'
                PVM_ObjOrderScheme = tag{2};
                if strcmp(PVM_ObjOrderScheme,'Interlaced')
                    pars.PVM_ObjOrderScheme = 3;
                elseif strcmp(PVM_ObjOrderScheme,'Sequential')
                    pars.PVM_ObjOrderScheme = 1;
                end
        end
    end
    fclose(fid);
 
    
    %Read reco 
    path_reco = deblank([pathscan filesep 'pdata' filesep '1' filesep 'reco']);
    fid       = fopen(path_reco,'rt');
    while feof(fid) == 0        
            line    =   fgetl(fid);
            tag     =	strread(line,'%s','delimiter','=');
            switch tag{1}
                case '##$RECO_maxima'
                    temp = fgetl(fid);
                    maxima = strread(temp,'%f');
                    pars.RECO_maxima = max(maxima);
                case '##$RECO_map_max'
                    temp = fgetl(fid);
                    maxima = strread(temp,'%f');
                    pars.RECO_map_maxima = max(maxima);    
                case '##$RECO_fov'
                    if strcmp(tag{2},'( 2 )')
                        [a, b]               =   strread(fgetl(fid),'%f %f');
                        pars.FOV(1)              =   10*a; 
                        pars.FOV(2)              =   10*b;                        
                    elseif strcmp(tag{2},'( 3 )')
                        [a, b, c]               =   strread(fgetl(fid),'%f %f %f');
                        pars.FOV(1)              =   10*a; 
                        pars.FOV(2)              =   10*b; 
                        pars.FOV(3)              =   10*c;
                    end
                case '##$RECO_size'
                    if strcmp(tag{2},'( 2 )')
                        [pars.dims(1) pars.dims(2)]   =   strread(fgetl(fid),'%d %d');    
                    elseif strcmp(tag{2},'( 3 )')
                        [pars.dims(1) pars.dims(2) pars.dims(3)]   =   strread(fgetl(fid),'%d %d %d');
                        pars.idist               =   pars.resol(3);
                    end                    
                case '##$RECO_wordtype'
                    pars.tp  =   tag{2};
                    switch pars.tp
                        case '_8BIT_USGN_INT'
                            pars.tp  =   'uint8';
                            pars.datatype = 2;
                            pars.bitpix   = 8;
                        case '_16BIT_SGN_INT' 
                            pars.tp  =   'int16'; 
                            pars.datatype = 4;
                            pars.bitpix   = 16;
                        case '_32BIT_SGN_INT' 
                            pars.tp  =   'int32'; 
                            pars.datatype = 8;
                            pars.bitpix   = 32;
                        case '_32BIT_FLT' 
                            pars.tp  =   'float32';
                            pars.datatype = 16;
                            pars.bitpix   = 32;
                        case '_64BIT_FLT' 
                            pars.tp  =   'float64'; 
                            pars.datatype = 64;
                            pars.bitpix   = 64;
                        case '_8BIT_SGN_INT' 
                            pars.tp  =   'int8'; 
                            pars.datatype = 256;
                            pars.bitpix   = 8;
                        case '_16BIT_USGN_INT' 
                            pars.tp  =   'uint16'; 
                            pars.datatype = 512;
                            pars.bitpix   = 16;
                        case '_32BIT_USGN_INT' 
                            pars.tp  =   'uint32';  
                            pars.datatype = 768;
                            pars.bitpix   = 32;            
                    end
                case '##$RECO_byte_order'
                    if strcmp (tag{2},'littleEndian')
                        pars.endian  =   'l';
                    elseif strcmp (tag{2},'bigEndian')
                        pars.endian  =   'b';                        
                    end
                case '##$RecoNumInputChan'
                    pars.n_coils     =   str2num(tag{2});
                case '##$RECO_image_type'
                    if strmatch(tag{2},'COMPLEX_IMAGE')
                        pars.cmpx    =   1;
                    end
                case '##$RecoScaleChan'
                    scales        =   [];
                    while true
                        read        =   fgetl(fid);
                        header      =   strread(read,'%c','delimiter','#');
                        if strcmp(header(1),'#') ||strcmp(header(1),'$') break; end
                        scales   =   [scales;strread(read,'%f','delimiter','\\ ')];
                    end
                    pars.scale =   scales;  
            end  
    end 
    fclose(fid);
    
    % Spatial resolution
    pars.resol(1:2)     =   [pars.FOV(1)/pars.dims(1); pars.FOV(2)/pars.dims(2)];
    if pars.acq_dim==2
        pars.FOV(3)     =   pars.dims(3)*pars.idist;
    end
    pars.dims           =   cast(pars.dims,'int16');    
    
    % Position calculation. FOV, resolution order
    dim             =   pars.dims;

    switch pars.orient
        case 'axial' 
            if strcmp(pars.r_out,'L_R')                     
                pars.vect   =   [1,2,3];              
            elseif   strcmp(pars.r_out,'A_P')
                pars.vect   =   [2,1,3];               
            end
        case 'sagittal'
            if strcmp(pars.r_out,'H_F')                 
                pars.vect   =   [2,1,3];            
            elseif strcmp(pars.r_out,'A_P')
                pars.vect   =   [1,2,3];               
            end
        case 'coronal'
            if strcmp(pars.r_out,'H_F')                 
                pars.vect   =   [2,1,3];             
            elseif strcmp(pars.r_out,'L_R')
                pars.vect      =   [1,2,3];               
            end        
    end
    pars.FOV    =   pars.FOV(pars.vect);                  
    half_vx     =   pars.resol(pars.vect)/2;  
    pars.offset =   [pars.offset(pars.vect(1));pars.offset(pars.vect(2));-pars.offset(pars.vect(3))]; 
    m_or2       =   [pars.m_or(:,pars.vect(1)) pars.m_or(:,pars.vect(2)) pars.m_or(:,pars.vect(3))]; 
                
    if isfield(pars,'MBfactor') 
        pars.FOV(3) = pars.FOV(3)*pars.MBfactor;
    end
    shift       =   -pars.FOV/2+half_vx-pars.offset;
    pars.pos0   =   m_or2*shift;
    pars.m_or   =   m_or2;
end


function im = ft2d(kdata)
    kdata = ifftshift(ifftshift(kdata,1),2);
    im    = fftshift(fftshift(fft2(kdata),1),2)/size(kdata,1)/size(kdata,2);
end

function im = ift2d(kdata)
    kdata = ifftshift(ifftshift(kdata,1),2);
    im    = fftshift(fftshift(ifft2(kdata),1),2);
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

% Slice-GRAPPA reconstruction
function [kspace_recon,sosimg] = sliceGRAPPA(kspace_und, kRef, kSize)

kRefCombined = sum(kRef,5);
if nargin < 3
    kSize  = [11 11];
end

lambda = 0.01;

kSize   = kSize + [1 1] - mod(kSize,2);
winSize = floor(kSize/2);

[Matx,Maty,Matz,nCoil,frames] = size(kspace_und);
kspace_recon = zeros(size(kspace_und));

for slice = 1:Matz
    kernel = calibrate_sliceGRAPPA(kRefCombined(:,:,slice,:),kRef(:,:,slice,:,1),kSize, nCoil, lambda);
    for v = 1:frames
        kData  = squeeze(kspace_und(:,:,slice,:,v));
        for x = winSize(1)+1:Matx-winSize(1)
            for y = winSize(2)+1:Maty-winSize(2)
                temp = kData(x-winSize(1):x+winSize(1),y-winSize(2):y+winSize(2),:);
                kspace_recon(x,y,slice,:,v) = reshape(temp(:).'*kernel,1,1,1,[]);
            end
        end
    end
end

sosimg = squeeze(sos(ft2d(kspace_recon),4));
end

% sum-of-squares
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The following functions were modified from the ESPIRiT toolbox by Michael Lustig
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calibrate slice-GRAPPA kernel
function [kernel,rawkernel] = calibrate_sliceGRAPPA(kRefCombined,kRef,kSize, nCoil, lambda)

    [sx,sy,nc] = size(kRefCombined);

    tmp1 = im2row(kRefCombined,kSize); [tsx,tsy,tsz] = size(tmp1);
    A   = reshape(tmp1,tsx,tsy*tsz);

    tmp2 = im2row(kRef,kSize); [tsx,tsy,tsz] = size(tmp2);
    y   = reshape(tmp2,tsx,tsy*tsz);

    dummyK = zeros(kSize(1),kSize(2),nCoil); dummyK(ceil((end+1)/2),ceil((end+1)/2),:) = 1;
    idxY = find(dummyK);
    y    = y(:,idxY);

    AtA = A'*A;
    Aty = A'*y;

    lambda    = norm(AtA,'fro')/size(AtA,1)*lambda;
    rawkernel = (AtA + eye(size(AtA))*lambda)\Aty;
    kernel    = Aty*0;
    kernel(:) = rawkernel; 

end

% Generate calibration rows
function res = im2row(im, winSize)

    [sx,sy,sz] = size(im);

    res = zeros((sx-winSize(1)+1)*(sy-winSize(2)+1),prod(winSize),sz);
    count=0;
    for y=1:winSize(2)
        for x=1:winSize(1)
            count = count+1;
            res(:,count,:) = reshape(im(x:sx-winSize(1)+x,y:sy-winSize(2)+y,:),...
                (sx-winSize(1)+1)*(sy-winSize(2)+1),1,sz);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

