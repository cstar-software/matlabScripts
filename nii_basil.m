function exitCode = nii_basil(asl, t1, aslRev, inCalScan, anatDir, dryRun, overwrite, outDir)

  root = '/Volumes/Ryan/aslRorden';
  addpath([root '/dicm2nii']);
  addpath([root '/spmScripts']);

  % TESTING
  TESTING = ~exist('asl','var');
  
  % handle overwritting
  % if the output directory contains a logfile then
  % assume the data has already been processed
  if exist('overwrite','var') && ~overwrite
    log = [outDir '/logfile'];
    if exist(log, 'file')
      disp([outDir ' is already processed']);
      exitCode = 0;
      return
    else
      disp(['resume processing ' outDir]);
    end
  end

  % v = version('-release');
  % disp(['MATLAB VERSION: ' v]);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % LINKS:
  %   https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/BASIL/UserGuide
  %   https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/BASIL/Tutorial
  %   https://users.fmrib.ox.ac.uk/~chappell/asl_primer/data/index.html (ASL Data for Examples)

  % pASL (Siemens): "PulseSequenceDetails" = "ep2d_pasl",
  %   https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=ind1404&L=FSL&D=0&P=185303
  %   set "labeling" to pASL
  %   Bolus duration: inversion time 1
  %   TIs: inversion time 2

  % pcASL (Oxford): "PulseSequenceDetails" = "ep2d_VEPCASL",
  %   set "labeling" to cASL/pcASL
  %   Bolus duration: BolusDuration
  %   PLDs: InitialPostLabelDelay
  %

  % pcASL (U Southern California/Danny Wang): "PulseSequenceDetails" = "ep2d_pcasl_UI_PHC",
  %   set "labeling" to cASL/pcASL
  %   Bolus duration: NumRFBlocks*0.0185*1000 %0.0185ms, convert to sec %https://www.mccauslandcenter.sc.edu/crnl/tools/asl
  %   PLD: PostLabelDelay
  %   No calibration image

  % pASL (FAIREST): "PulseSequenceDetails": "ep2d_fairest_UI_iPAT"
  %   http://www.pubmed.com/21606572, 
  %   set "labeling" to cASL/pcASL
  %   Bolus duration: PostInversionDelay
  %   PLD: PostLabelDelay
  %   Creates 2 series: M0 Calibration Scan (Proton Density) and label+tagged pairs
  % 
  %   TODO: SliceTiming: there is no slice timing! (use from another
  %   sequence)

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  % GLOBALS

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
  SEQUENCE_PASL_SIEMENS = 'ep2d_pasl';
  SEQUENCE_PCASL_OXFORD = 'ep2d_VEPCASL';
  SEQUENCE_PCASL_USC_WANG = 'ep2d_pcasl_UI_PHC';
  SEQUENCE_PASL_FAIREST = 'ep2d_fairest_UI_iPAT';
  
  % TODO:
  % - are these hard coded values correct for all sequences?
      % '--t1b 1.65 ',           ... 
      % '--alpha 0.85 ',         ...   

  % EXAMPLE SEQUNCES FOR TESTING
  if TESTING
    disp ' '; 
    disp '========= TESTING =========';
    disp ' ';
    
    % chris (SEQUENCE_PCASL_OXFORD)
    % sequence = SEQUENCE_PCASL_OXFORD;
    % outDir = '/Volumes/Ryan/aslRorden/out';
    % dir = '/Volumes/Ryan/aslRorden';
    % asl = [dir '/24_asl_ep2d_PCASL.nii'];
    % t1 = [dir '/6_anat_T1w.nii'];

    % william (SEQUENCE_PCASL_OXFORD)
    % outDir = '/Volumes/Ryan/aslRorden/PCASL_20190208135617';
    % dir = outDir;
    % asl = [dir '/PCASL_20190208135617_to_ep2d_PCASL_20190208135617_11.nii.gz'];
    % aslRev = [dir '/PCASL_20190208135617_to_ep2d_PCASL_PA_20190208135617_15.nii.gz'];
    % t1 = [dir '/PCASL_20190208135617_T1_mprage_ns_sag_p2_iso_1mm_192_20190208135617_2.nii.gz'];

    % SEQUENCE_PCASL_USC_WANG
    % outDir = '/Volumes/Ryan/aslRorden/pasl_pcasl_old/pcasl/20120607135306';
    % dir = outDir;
    % asl = [dir '/20120607135306_pCASL_1200ms_3.5sTR_10.nii.gz'];
    % t1 = [dir '/20120607135306_T1_mprage_ns_sag_p2_iso_1.0mm_192_7.nii.gz'];

    % SEQUENCE_PASL_SIEMENS
    % outDir = '/Volumes/Ryan/aslRorden/pasl_pcasl_old/pasl';
    % dir = outDir;
    % asl = [dir '/4_ep2d_tra_pasl_p2.nii'];
    % t1 = [dir '/2_Anatomical_1x1x1_ipat=2_nw.nii'];

    % SEQUENCE_PASL_FAIREST
    % outDir = '/Volumes/Ryan/aslRorden/pasl_fairest_new';
    % dir = outDir;
    % asl = [dir '/20091031130655_PASL_fairest_16x5_no-fs1.5_7.nii.gz'];
    % t1 = [dir '/20091031130655_Anatomical_1x1x1_ipat=2_nw_4.nii.gz'];
    % inCalScan = [dir '/20091031130655_PASL_fairest_16x5_M0_nofs_5.nii.gz'];
    
    % sen pcasl
    % outDir = '/Volumes/Ryan/aslRorden/sen';
    % dir = outDir;
    % asl = [dir '/16_ep2d_pcasl_UI_PHC_p4.nii'];
    % t1 = [dir '/4_T1_MPRage.nii'];

    outDir = '/Volumes/Ryan/aslRorden/sen';
    dir = '/Volumes/Chris5TB/SEN/nii/2035';
    asl = [dir '/Sen2039_20141024090739_4_ep2d_pcasl_UI_PHC_p4.nii'];
    t1 = [dir '/Sen2039_20141024090739_3_T1_MPRage.nii'];

    aslOutput = [outDir '/oxford_asl_output'];
    anatDir = '';
    dryRun = true;
  else
    % TODO: is this ok for m0 files that are made?
    aslOutput = outDir;
  end

  if dryRun
    disp '========= DRY RUN =========';
    disp ' ';
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



  % check if input files exist or prompt to find new ones
  % if ~exist('aslRev','var') && ~exist('asl','var')
  %    [A,Apth] = uigetfile({'*.nii;*.gz;';'*.*'},'Select reversed phase ASL image');
  %    aslRev = [Apth, A];
  % elseif ~exist('aslRev','var')
  %   aslRev = [];
  % end
  % if ~exist('asl','var')
  %    [A,Apth] = uigetfile({'*.nii;*.gz;';'*.*'},'Select ASL image');
  %    asl = [Apth, A]; 
  % end
  % if ~exist('t1','var')
  %    [A,Apth] = uigetfile({'*.nii;*.gz;';'*.*'},'Select T1 image');
  %    t1 = [Apth, A];
  % end

  if ~exist('aslRev','var')
    aslRev = [];
  end

  % verify inputs
  if ~exist(t1,'file'), error('Requires T1 scan %s\n', t1); end
  if ~exist(asl,'file'), error('Requires %s\n', asl); end
  if ~isempty(aslRev) && ~exist(aslRev,'file'), error('Requires %s\n', aslRev); end
  if ~exist(asl,'file'), error('Requires %s\n', asl); end

  if ~exist('outDir','var') outDir = pwd; end
  if ~exist(outDir,'dir') mkdir(outDir); end
  if ~exist(aslOutput,'dir') mkdir(aslOutput); end

  % will create:
  aslX = [];        %asl data with M0 removed
  aslRevX = [];     %reverse phased first volulme
  m0X = [];         %M0 calibration scan
  if exist(aslRev,'file')
      nii = nii_tool('load', aslRev);
      [~,n,x] = nii_fileparts(aslRev);
      niiRev.hdr = nii.hdr;
      niiRev.hdr.dim(1) = 3; %3d not 4D data
      niiRev.hdr.dim(5) = 1; %single volume
      aslRevX = fullfile(outDir, ['m0', n, x]);
      niiRev.img = nii.img(:,:,:,1);
      nii_tool('save', niiRev, aslRevX);
  end
  
  % load ASL
  if ~exist(asl,'file'), error('Unable to find %s\n', asl); end
  nii = nii_tool('load', asl);
  [TR, PLDs, BD, SecPerSlice, EffectiveEchoSpacing, sequence] = readJson(asl);

  % check number of volumes
  vols = nii.hdr.dim(5);
  if mod( (vols-1), numel(PLDs)), error('%d volumes and %d PLDs, expected 1 calibration image with remaining volumes divisble by PLDs', vols, numel(PLDs)); end
  repeats = (vols-1) / numel(PLDs);
  if (repeats < 2) , error('%d volumes and %d PLDs, expected multiple volumes per PLD', vols, numel(PLDs)); end
  if (repeats < 2) , error('%d volumes and %d PLDs, expected multiple volumes per PLD', vols, numel(PLDs)); end
  
  % do some additional formatting depending on sequence type
  switch sequence
    case {SEQUENCE_PASL_SIEMENS, SEQUENCE_PCASL_OXFORD}
      if (mod(repeats, 2) ~= 0), error('%d volumes, %d repeats and %d PLDs, expected even number of volumes per PLD (same number label and control)', vols, repeats, numel(PLDs)); end 
      %save calibration image
      [~,n,x] = nii_fileparts(asl);
      niiM0.hdr = nii.hdr;
      niiM0.hdr.dim(1) = 3; %3d not 4D data
      niiM0.hdr.dim(5) = 1; %single volume
      m0X = fullfile(outDir, ['m0_', n, x]);
      niiM0.img = nii.img(:,:,:,1);
      nii_tool('save', niiM0, m0X);
      %save subsequent volumes (ASL images)
      [~,n,x] = nii_fileparts(asl);
      niiAsl.hdr = nii.hdr;
      niiAsl.hdr.dim(5) = nii.hdr.dim(5) - 1; %single volume
      aslX = fullfile(outDir, ['asl_', n, x]);
      niiAsl.img = nii.img(:,:,:,2:end);
      nii_tool('save', niiAsl, aslX);
    case SEQUENCE_PCASL_USC_WANG
      aslX = asl;
    case SEQUENCE_PASL_FAIREST
      aslX = asl;
      if ~exist(inCalScan, 'file')
        error('SEQUENCE_PASL_FAIREST requires a calibration image to be passed in directly');
      end
      m0X = inCalScan;
    otherwise
      error('invalid sequence %d', sequence)
  end

  % Grouping order:
  %   repeats: ibf=tis
  %   label/control pairs: ibf=rpt

  % Label/Control pairs: --iaf
  %   label then control: (iaf=tc)
  %   control then label: (iaf=ct)

  % sequence specific options
  switch sequence
    case SEQUENCE_PCASL_OXFORD
      options.cASL = true;
      options.groupingOrder = 'rpt';
      options.labelControlPairs = 'tc';
    case SEQUENCE_PASL_SIEMENS
      options.cASL = false;
      options.groupingOrder = 'rpt';
      options.labelControlPairs = 'tc';
    case SEQUENCE_PCASL_USC_WANG
      options.cASL = true;
      options.groupingOrder = 'rpt';
      options.labelControlPairs = 'tc';
    case SEQUENCE_PASL_FAIREST
      options.cASL = false;
      % TODO: what grouping order is FAIREST? chris seems to think
      % label/control pairs
      options.groupingOrder = 'tis';
      options.labelControlPairs = 'tc';
  end

  options.output = aslOutput;
  options.inputImage = aslX;
  options.calibrationImage = m0X;
  options.reversedCalibrationImage = aslRevX;
  options.T1Image = t1;
  options.TR = TR;
  options.PLDs = PLDs;
  options.BD = BD;
  options.secPerSlice = SecPerSlice;
  options.effectiveEchoSpacing = EffectiveEchoSpacing;
  
  % generate anat .struct file
  if isempty(anatDir)
    anatDir = options.output;
  end
  options.anatFile = [anatDir '/struc.anat'];

  if ~exist(options.anatFile, 'file')
    command = ['fsl_anat -i "' options.T1Image '" -o ' anatDir '/struc'];
    if ~dryRun
      exitCode = fslCmd(command);
      if exitCode ~= 0
        return;
      end
    else
      % show command for dry runs
      disp(command);    
    end
  elseif dryRun
    disp (['found .anat file at "' options.anatFile '"']);
  end

  command = makeCommand(options);

  % run the command
  if ~isempty(command) && ~dryRun
    exitCode = fslCmd(command);
  else
    exitCode = 0;
  end

  % ============================================================================================================

  function command = makeCommand(options)

    % generate tis
    tis = '';
    for n = 1:numel(options.PLDs)
      % TODO: is this true for all sequences?
      label = options.BD + options.PLDs(n);
      tis = strcat(tis, num2str(label)); 
      if n < 6
          tis = strcat(tis, ',');
      end
    end

    % program
    part1 = 'oxford_asl';

    % input data
    % tis = PostLabelDelay+BolusTime
    % casl = labeling cASL/pcASL
    if options.cASL
      cASLFlag = '--casl';
    else
      cASLFlag = '';
    end

    part2 = [
      '-i "', options.inputImage, '" ',             ...  % input image
      cASLFlag, ' ',                                ...  % Labeling (pASL or cASL/pcASL)
      '--ibf=', options.groupingOrder, ' ',         ...  % (Grouping order): Input block format. Specifically for multi-delay (multi-PLD) ASL data to identify whther the individual delays/PLDs are groups togther or by repeats of the same sequence of PLDs.
      '--iaf=', options.labelControlPairs, ' ',     ...  % (Label/Control pairs): Input ASL format: specifies if the data has already been label-control subtracted (diff, default), or is in the form of label(tag)-control pairs (tc or ct depending upon if label/tag is first).
      '--tis ', tis, ' ',                           ...  % TIs/PLDs depending on labeling (pASL/pcASL)
      '--bolus ', num2str(options.BD), ' ',         ...  % Bolus duration
      '--slicedt ', num2str(options.secPerSlice)         % Time per slice (ms)
    ];

    % structural
    if isempty(options.anatFile)
      part3 = []; %no structural file
    else
      part3 = ['--fslanat="', options.anatFile, '"'];
    end

    % calibration
    if isempty(options.calibrationImage)
      part4 = []; %no calibration image
    else
      part4 = [
        '-c "', options.calibrationImage, '" ',   ...
        '--tr ', num2str(options.TR), ' ',        ...
        '--cgain 1.00 ',                          ...
        '--cmethod voxel'
      ];
    end

    % distortion correction
    if isempty(options.reversedCalibrationImage)
      part5 = []; %no phase reversed image
    else
      part5 = [
        '--cblip="',  options.reversedCalibrationImage ,'" ',                   ...
        '--echospacing=', num2str(options.effectiveEchoSpacing * 1000) ,' '     ...
        '--pedir=y'
      ];
    end

    % analysis
    part6 = [
      '--wp ',                 ...  % analysis conforms to white paper
      '--t1b 1.65 ',           ...  
      '--alpha 0.85 ',         ...  % inversion efficiency
      '--fixbolus ',           ...  % fixed label duration
      '--spatial ',            ...  % adaptive spatial regularization
      '--artoff ',             ...
      '--pvcorr ',             ...
      '--mc ',                 ...  % motion correction
      '-o "', options.output, '"'
    ];

    command = [part1 ' ' part2 ' ' part3 ' ' part4 ' ' part5 ' ' part6]
  end

  function [TR, PLDs, BD, SecPerSlice, EffectiveEchoSpacing, sequence] = readJson(niiName)
    [p,n,x] = nii_fileparts(niiName);
    json = fullfile(p,[n,'.json']);
    if ~exist(json,'file'), error('Unable to find %s\n', json); end
    js = loadJson(json);
    if isfield(js, 'MultibandAccelerationFactor'), error('Update to support multiband'); end
    if ~isfield(js, 'PulseSequenceDetails'), error('Require PulseSequenceDetails'); end

    % get sequence name
    result = regexp(js.PulseSequenceDetails, '^%CustomerSeq%_(\w+)$', 'tokens');
    sequence = result{1}{1};

    switch sequence

      case SEQUENCE_PCASL_OXFORD
        if ~isfield(js, 'InversionTime'), error('Require InversionTime'); end
        if ~isfield(js, 'BolusDuration'), error('Require BolusDuration'); end
        if ~isfield(js, 'InitialPostLabelDelay'), error('Require InitialPostLabelDelay'); end
        if ~isfield(js, 'RepetitionTime'), error('Require RepetitionTime'); end
        if ~isfield(js, 'SliceTiming'), error('Require SliceTiming'); end
        if ~isfield(js, 'EffectiveEchoSpacing'), error('Require EffectiveEchoSpacing'); end
        if iscell(js.SliceTiming)
          SecPerSlice = (max(cell2mat(js.SliceTiming(:))) -min(cell2mat(js.SliceTiming(:))))/(numel(js.SliceTiming)-1);
        else
          SecPerSlice = (max(js.SliceTiming(:)) -min(js.SliceTiming(:)))/(numel(js.SliceTiming)-1);
        end
        % TI is not used anywhere so we removed it
        % TI = js.InversionTime;
        BD = js.BolusDuration;
        PLDs = js.InitialPostLabelDelay;
        TR = js.RepetitionTime;
        EffectiveEchoSpacing = js.EffectiveEchoSpacing;

      case SEQUENCE_PASL_SIEMENS
        if ~isfield(js, 'InversionTime1'), error('Require InversionTime1'); end
        if ~isfield(js, 'InversionTime2'), error('Require InversionTime2'); end
        if ~isfield(js, 'RepetitionTime'), error('Require RepetitionTime'); end
        if ~isfield(js, 'SliceTiming'), error('Require SliceTiming'); end
        if ~isfield(js, 'EffectiveEchoSpacing'), error('Require EffectiveEchoSpacing'); end
        if iscell(js.SliceTiming)
          SecPerSlice = (max(cell2mat(js.SliceTiming(:))) -min(cell2mat(js.SliceTiming(:))))/(numel(js.SliceTiming)-1);
        else
          SecPerSlice = (max(js.SliceTiming(:)) -min(js.SliceTiming(:)))/(numel(js.SliceTiming)-1);
        end

        BD = js.InversionTime1 / 1000;
        PLDs = [js.InversionTime2 / 1000];
        TR = js.RepetitionTime;
        EffectiveEchoSpacing = js.EffectiveEchoSpacing;

      case SEQUENCE_PCASL_USC_WANG
        if ~isfield(js, 'NumRFBlocks'), error('Require NumRFBlocks'); end
        if ~isfield(js, 'PostLabelDelay'), error('Require PostLabelDelay'); end
        if ~isfield(js, 'RepetitionTime'), error('Require RepetitionTime'); end
        if ~isfield(js, 'SliceTiming'), error('Require SliceTiming'); end
        if ~isfield(js, 'EffectiveEchoSpacing'), error('Require EffectiveEchoSpacing'); end
        if iscell(js.SliceTiming)
          SecPerSlice = (max(cell2mat(js.SliceTiming(:))) -min(cell2mat(js.SliceTiming(:))))/(numel(js.SliceTiming)-1);
        else
          SecPerSlice = (max(js.SliceTiming(:)) -min(js.SliceTiming(:)))/(numel(js.SliceTiming)-1);
        end

        BD = js.NumRFBlocks * 0.0185;
        PLDs = [js.PostLabelDelay];
        TR = js.RepetitionTime;
        EffectiveEchoSpacing = js.EffectiveEchoSpacing;

      case sSEQUENCE_PASL_FAIREST
        if ~isfield(js, 'PostInversionDelay'), error('Require PostInversionDelay'); end
        if ~isfield(js, 'PostLabelDelay'), error('Require PostLabelDelay'); end
        if ~isfield(js, 'RepetitionTime'), error('Require RepetitionTime'); end
        if ~isfield(js, 'EffectiveEchoSpacing'), error('Require EffectiveEchoSpacing'); end
        %   TODO: SliceTiming is not available, hardcode something meaningful
        SecPerSlice = 1;
        BD = js.PostInversionDelay / 1000;
        PLDs = [js.PostLabelDelay / 1000];
        TR = js.RepetitionTime;
        EffectiveEchoSpacing = js.EffectiveEchoSpacing;

      otherwise
        error('Invalid sequence %d', sequence);
    end
  end

  function [json, fnm] = loadJson(fnm)
    p = fileparts(which(mfilename));
    if ~exist('fnm','var'), fnm = fullfile(p, 'config.json'); end
    if ~exist(fnm,'file')
      fprintf('Creating empty JSON: Unable to find %s\n', fnm);
      json = [];
      return;
    end
    fid = fopen(fnm);
    raw = fread(fid,inf);
    str = char(raw');
    fclose(fid);
    if exist('jsondecode', 'builtin')
      json = jsondecode(str);
    else
      % https://www.mathworks.com/matlabcentral/fileexchange/20565-json-parser
      js = parse_json(str);
      json = js{1};
    end
  end

  function path = findBestPath(paths)
    path = nan;
    for i = 1:numel(paths)
      if exist(paths{i},'dir')
        path =  paths{i};
        return
      end
    end
  end

  function exitCode = fslCmd(Cmd)
    % find fsl path
    fsldir = findBestPath({ '/usr/local/fsl/',
                            '/Volumes/Ryan/fsl/',
                            '/usr/share/fsl/6.0/'
                            });
    if ~exist(fsldir,'dir')
      error('%s: fsldir not found',mfilename);
    end
    setenv('FSLDIR', fsldir);
    flirt = [fsldir 'bin/flirt'];
    if ~exist(flirt,'file')
      error('%s: flirt (%s) not found',mfilename,flirt);
    end
    command=sprintf('sh -c ". %setc/fslconf/fsl.sh; %sbin/%s"\n',fsldir,fsldir, Cmd);
    %command=sprintf('sh -c ". %setc/fslconf/fsl.sh; %sbin/%s" >log.txt\n',fsldir,fsldir, Cmd);

    fprintf(command);
    exitCode = system(command);
  end

end