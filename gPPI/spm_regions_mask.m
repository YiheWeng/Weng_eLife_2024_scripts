function [Y,xY] = spm_regions_mask(job)
% VOI time-series extraction of adjusted data (& local eigenimage analysis)
% FORMAT [Y,xY] = spm_regions(xSPM,SPM,hReg,[xY])
%
% xSPM   - structure containing specific SPM, distribution & filtering details
% SPM    - structure containing generic analysis details
% hReg   - Handle of results section XYZ registry (see spm_results_ui.m)
%
% Y      - first scaled eigenvariate of VOI {i.e. weighted mean}
% xY     - VOI structure
%       xY.xyz          - centre of VOI {mm}
%       xY.name         - name of VOI
%       xY.Ic           - contrast used to adjust data (0 - no adjustment)
%       xY.Sess         - session index
%       xY.def          - VOI definition
%       xY.spec         - VOI definition parameters
%       xY.str          - VOI description as a string
%       xY.XYZmm        - Co-ordinates of VOI voxels {mm}
%       xY.y            - [whitened and filtered] voxel-wise data
%       xY.u            - first eigenvariate {scaled - c.f. mean response}
%       xY.v            - first eigenimage
%       xY.s            - eigenvalues
%       xY.X0           - [whitened] confounds (including drift terms)
%
% Y and xY are also saved in VOI_*.mat in the SPM working directory.
% (See spm_getSPM for details on the SPM & xSPM structures)
%
% FORMAT [Y,xY] = spm_regions('Display',[xY])
%
% xY     - VOI structure or filename
%
%__________________________________________________________________________
%
% spm_regions extracts a representative time course from voxel data in
% terms of the first eigenvariate of the filtered and adjusted response in
% all suprathreshold voxels within a specified VOI centred on the current
% MIP cursor location. Responses are adjusted by removing variance that
% can be predicted by the null space of the F contrast specified (usually 
% an F-contrast testing for all effects of interest).
%
% If temporal filtering has been specified, then the data will be filtered.
% Similarly for whitening. Adjustment is with respect to the null space of
% a selected contrast, or can be omitted.
%
% For a VOI of radius 0, the [adjusted] voxel time-series is returned, and
% scaled to have a 2-norm of 1. The actual [adjusted] voxel time series can
% be extracted from xY.y, and will be the same as the [adjusted] data 
% returned by the plotting routine (spm_graph.m) for the same contrast.
%__________________________________________________________________________
% Copyright (C) 1999-2015 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_regions.m 6436 2015-05-14 10:05:27Z guillaume $

%--------------------------------------------------------------------------
load(cell2mat(job.spmmat));

%-Initialise VOI voxels coordinates
%--------------------------------------------------------------------------
[x,y,z] = ndgrid(1:SPM.xVol.DIM(1),1:SPM.xVol.DIM(2),1:SPM.xVol.DIM(3));  % 3D dim
XYZ     = [x(:),y(:),z(:)]'; clear x y z   % native space
XYZmm   = SPM.xVol.M(1:3,:) * [XYZ;ones(1,size(XYZ,2))];  % MNI coordinater

%-Estimate VOIs
%--------------------------------------------------------------------------
voi     = cell(1,numel(job.roi));
for i=1:numel(job.roi)
    voi = roi_estim(XYZmm,i,job,SPM,voi);
end

%-Evaluate resulting VOI
%--------------------------------------------------------------------------
voi     = roi_eval(voi,job.expression);  % mask information,job.expression: i1+i2

xSPM.XYZmm = XYZmm;
xSPM.XYZ   = XYZ;
xSPM.M     = SPM.xVol.M; % irrelevant here
 
% %-Get adjustment options and VOI name
% %--------------------------------------------------------------------------
%-Extract VOI time-series
%--------------------------------------------------------------------------
xY.name    = job.name;
xY.Ic      = job.adjust;
xY.Sess    = job.session;
xY.xyz     = []'; % irrelevant here
xY.def     = 'mask';
% xY.spec    = Vm;

if ~isfield(xY,'Ic')
    q(1)   = 0;
    Con    = {'<don''t adjust>'};
    q(2)   = NaN;
    Con{2} = '<adjust for everything>';
    for i = 1:length(SPM.xCon)
        if strcmp(SPM.xCon(i).STAT,'F')
            q(end + 1) = i;
            Con{end + 1} = SPM.xCon(i).name;
        end
    end
    if numel(Con) == 2
        warning('No F-contrast has been defined: are you sure?');
    end
    i     = spm_input('adjust data for (select contrast)','!+1','m',Con);
    xY.Ic = q(i);
end
 
%-If fMRI data then ask user to select session
%--------------------------------------------------------------------------
if isfield(SPM,'Sess') && ~isfield(xY,'Sess')
    s       = length(SPM.Sess);
    if s > 1
        s   = spm_input('which session','!+1','n1',s,Inf);
    end
    xY.Sess = s;
end

%-Specify VOI
%--------------------------------------------------------------------------
Q          = (find(voi>0))';
xY.XYZmm   = xSPM.XYZmm(:,Q);  % Mask coordinate

if isempty(xY.XYZmm)
    warning('Empty region.');
    Y = [];
    return;
end

%-Perform time-series extraction to all sessions if Inf is entered
%--------------------------------------------------------------------------
if isfield(SPM,'Sess') && isfield(xY,'Sess') && isinf(xY.Sess)
    if length(SPM.Sess) == 1
        xY.Sess = 1;
    else
        for i=1:length(SPM.Sess)
            xY.Sess = i;
            [tY{i},txY(i)] = spm_regions(xSPM,SPM,hReg,xY);
        end
        Y = tY; xY = txY;
        return;
    end
end
 
%-Extract required data from results files
%==========================================================================
spm('Pointer','Watch')
 
%-Get raw data, whiten and filter 
%--------------------------------------------------------------------------
y        = spm_data_read(SPM.xY.VY,'xyz',xSPM.XYZ(:,Q));  % mask TC
y        = spm_filter(SPM.xX.K,SPM.xX.W*y);               % Filter TC
 
 
%-Computation
%==========================================================================
 
%-Remove null space of contrast
%--------------------------------------------------------------------------
if xY.Ic ~= 0
 
    %-Parameter estimates: beta = xX.pKX*xX.K*y
    %----------------------------------------------------------------------
    beta  = spm_data_read(SPM.Vbeta,'xyz',xSPM.XYZ(:,Q));
 
    %-subtract Y0 = XO*beta,  Y = Yc + Y0 + e
    %----------------------------------------------------------------------
    if ~isnan(xY.Ic)
        y = y - spm_FcUtil('Y0',SPM.xCon(xY.Ic),SPM.xX.xKXs,beta);
    else
        y = y - SPM.xX.xKXs.X * beta;
    end
 
end
 
%-Confounds
%--------------------------------------------------------------------------
xY.X0     = SPM.xX.xKXs.X(:,[SPM.xX.iB SPM.xX.iG]);  % constant?
 
%-Extract session-specific rows from data and confounds
%--------------------------------------------------------------------------
try
    i     = SPM.Sess(xY.Sess).row;
    y     = y(i,:);
    xY.X0 = xY.X0(i,:);
end
 
% and add session-specific filter confounds
%--------------------------------------------------------------------------
try
    xY.X0 = [xY.X0 SPM.xX.K(xY.Sess).X0]; % cosines (high-pass filter)
end
try
    xY.X0 = [xY.X0 SPM.xX.K(xY.Sess).KH]; % Compatibility check
end
 
%-Remove null space of X0
%--------------------------------------------------------------------------
xY.X0     = xY.X0(:,any(xY.X0));
 
 
%-Compute regional response in terms of first eigenvariate
%--------------------------------------------------------------------------
[m,n]   = size(y);
if m > n
    [v,s,v] = svd(y'*y);
    s       = diag(s);
    v       = v(:,1);
    u       = y*v/sqrt(s(1));
else
    [u,s,u] = svd(y*y');
    s       = diag(s);
    u       = u(:,1);
    v       = y'*u/sqrt(s(1));
end
d       = sign(sum(v));
u       = u*d;
v       = v*d;
Y       = u*sqrt(s(1)/n);
 
%-Set in structure
%--------------------------------------------------------------------------
xY.y    = y;
xY.u    = Y;
xY.v    = v;
xY.s    = s;
 
 
% %-Save
% %==========================================================================
% str = ['VOI_' xY.name '.mat'];
% if isfield(xY,'Sess') && isfield(SPM,'Sess')
%     str = sprintf('VOI_%s_%i.mat',xY.name,xY.Sess);
% end
% save(fullfile(SPM.swd,str),'Y','xY', spm_get_defaults('mat.format'))


% cmd = 'spm_regions(''display'',''%s'')';
% fprintf('   VOI saved as %s\n',spm_file(fullfile(SPM.swd,str),'link',cmd));
 
%-Reset title
%--------------------------------------------------------------------------
% set(Finter,'Name',header);
% spm('Pointer','Arrow')


%==========================================================================
% function voi = roi_estim(xyz,n,job,SPM,voi)
% voi = roi_eval(voi,expr)
% function [SPM, xSPM] = getSPM(s)
% function c = get_centre(xyz,n,job,SPM,voi)
%==========================================================================
function voi = roi_estim(xyz,n,job,SPM,voi)

if ~isempty(voi{n}), return; end

voi{n} = false(SPM.xVol.DIM');

Q      = ones(1,size(xyz,2));

switch char(fieldnames(job.roi{n}))
    
    case 'sphere'
    %----------------------------------------------------------------------
        c              = get_centre(xyz,n,job,SPM,voi);
        r              = job.roi{n}.sphere.radius;
        voi{n}(sum((xyz - c*Q).^2) <= r^2) = true;
    
    case 'box'
    %----------------------------------------------------------------------
        c              = get_centre(xyz,n,job,SPM,voi);
        d              = job.roi{n}.box.dim(:);
        voi{n}(all(abs(xyz - c*Q) <= d*Q/2)) = true;
    
    case 'spm'
    %----------------------------------------------------------------------
        if isempty(job.roi{n}.spm.spmmat{1})
            job.roi{n}.spm.spmmat = job.spmmat;
        end
        [SPM1,xSPM]    = getSPM(job.roi{n}.spm);
        voi1           = zeros(SPM1.xVol.DIM');
        voi1(sub2ind(SPM1.xVol.DIM',xSPM.XYZ(1,:),xSPM.XYZ(2,:),xSPM.XYZ(3,:))) = 1;
        XYZ            = SPM1.xVol.iM(1:3,:)*[xyz; Q];
        voi{n}(spm_sample_vol(voi1, XYZ(1,:), XYZ(2,:), XYZ(3,:),0) > 0) = true;
    
    case 'mask'
    %----------------------------------------------------------------------
        v              = spm_vol(job.roi{n}.mask.image{1});
        t              = job.roi{n}.mask.threshold;
        iM             = inv(v.mat);
        XYZ            = iM(1:3,:)*[xyz; Q];
        voi{n}(spm_sample_vol(v, XYZ(1,:), XYZ(2,:), XYZ(3,:),0) > t) = true;
    
    case 'label'
    %----------------------------------------------------------------------
        v              = spm_vol(job.roi{n}.label.image{1});
        l              = job.roi{n}.label.list;
        iM             = inv(v.mat);
        XYZ            = iM(1:3,:)*[xyz; Q];
        voi{n}(ismember(spm_sample_vol(v, XYZ(1,:), XYZ(2,:), XYZ(3,:),0),l)) = true;
    
end

%==========================================================================
function voi = roi_eval(voi,expr)
for i=1:numel(voi)
    eval(sprintf('i%d=voi{%d};',i,i));
end
try
    eval(['voi=' expr ';']);
catch
    error('The expression cannot be evaluated.');
end

%==========================================================================
function [SPM, xSPM] = getSPM(s)
xSPM.swd       = spm_file(s.spmmat{1},'fpath');
xSPM.Ic        = s.contrast;
xSPM.n         = s.conjunction;
xSPM.u         = s.thresh;
xSPM.thresDesc = s.threshdesc;
xSPM.k         = s.extent;
xSPM.title     = '';
xSPM.Im        = [];
if ~isempty(s.mask)
    xSPM.Im    = s.mask.contrast;
    xSPM.pm    = s.mask.thresh;
    xSPM.Ex    = s.mask.mtype;
end
[SPM,xSPM]     = spm_getSPM(xSPM);

%==========================================================================
function c = get_centre(xyz,n,job,SPM,voi)
t             = char(fieldnames(job.roi{n}));
c             = job.roi{n}.(t).centre(:);
mv            = char(fieldnames(job.roi{n}.(t).move));
if strcmp(mv,'fixed'), return; end

m             = job.roi{n}.(t).move.(mv).spm;
e             = job.roi{n}.(t).move.(mv).mask;
k             = union(roi_expr(e), m);
for i=1:numel(k)
    voi       = roi_estim(xyz,k(i),job,SPM,voi);
end
try
    if isempty(job.roi{m}.spm.spmmat{1})
        job.roi{m}.spm.spmmat = job.spmmat;
    end
catch
    error('The SPM index does not correspond to a Thresholded SPM ROI.');
end
[mySPM, xSPM] = getSPM(job.roi{m}.spm);
XYZmm         = xSPM.XYZmm;
XYZ           = SPM.xVol.iM(1:3,:)*[XYZmm;ones(1,size(XYZmm,2))];
Z             = xSPM.Z;
if ~isempty(e)
    R         = spm_sample_vol(uint8(roi_eval(voi,e)), ...
                  XYZ(1,:), XYZ(2,:), XYZ(3,:),0) > 0;
    XYZ       = XYZ(:,R);
    XYZmm     = xSPM.XYZmm(:,R);
    Z         = xSPM.Z(R);
end
[N, Z, M]     = spm_max(Z,XYZ);
if isempty(Z)
    warning('No voxel survived. Default to user-specified centre.');
    return
end

str           = '[%3.0f %3.0f %3.0f]';
switch mv
    case 'global'
        [i,j] = max(Z);
        nc    = SPM.xVol.M(1:3,:)*[M(:,j);1];
        str   = sprintf(['centre moved to global maximum ' str],nc);
    case 'local'
        XYZmm = SPM.xVol.M(1:3,:)*[M;ones(1,size(M,2))];
        nc    = spm_XYZreg('NearestXYZ',c,XYZmm);
        str   = sprintf(['centre moved from ' str ' to ' str],c,nc);
    case 'supra'
        nc    = spm_XYZreg('NearestXYZ',c,XYZmm);
        str   = sprintf(['centre moved from ' str ' to ' str],c,nc);
    otherwise
        error('Unknown option: ''%s''.',mv);
end
c             = nc;
fprintf(['   ' upper(t(1)) t(2:end) ' ' str '\n']);                     %-#

