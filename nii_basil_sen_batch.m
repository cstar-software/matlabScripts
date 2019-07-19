function exitCode = nii_basil_sen_batch(dataFolder)

  % constants
  ASL_NAME = 'ep2d_pcasl_UI_PHC_p4';
  T1_NAME = 'T1_MPRage';
  
  % TESTING  
  if ismac
    dataFolder = '/Volumes/Chris5TB/SEN/nii';  
  elseif isunix
    dataFolder = '/media/research/Chris5TB/SEN/nii';
  end
  
  if ~exist('dataFolder','var')
    [dataFolder] = uigetdir();
  end

  if ~exist(dataFolder,'dir')
    error('data folder "%s" not found', dataFolder);
  end

  outputFolder = [dataFolder '/output'];
  if ~exist(outputFolder,'dir')
    error('output folder "%s" not found', outputFolder);
  end

  % show some stuff
  disp(dataFolder);
  disp(outputFolder);
  
  total = 0;
  files = dir([dataFolder '/*'])
  for file = files'
    if file.name(1) ~= '.'
      path = [dataFolder '/' file.name];
      if isdir(path) && ~strcmp(path, outputFolder)
        exitCode = process_folder(path);
        if exitCode ~= 0
          return
        end
        total = total + 1;
      end
    end
  end
  
  % finished
  disp(['finished processing ' num2str(total) ' files.']);
  exitCode = 0;
  
  % ============================================================================================================

  function exitCode = process_folder(folder)
    files = dir([folder '/*']);
    group = [];
    currentID = '';
    count = 0;
    lastFile = false;
    for file = files'
      count = count + 1;
      if file.name(1) ~= '.'
        result = regexp(file.name, '^Sen(\w+)_(\d+)_(\d+)_(.*)\.(\w+)$', 'tokens');
        
        % invalid file name
        if isempty(result)
          error(['file name ', file.name, ' doesn''t match pattern, in ', folder]);
        end
        
        ext = char(result{1}{5});

        % ignore .json files
        if isequal(ext, 'json')
          continue
        end

        parts.path = [folder '/' file.name];
        parts.id = char(result{1}{1});
        parts.series = str2num(result{1}{3});
        parts.name = char(result{1}{4});
        parts.ext = ext;
        parts.isT1 = false;
        parts.isImage = false;
        if isequal(parts.name, T1_NAME)
          parts.isT1 = true;
        elseif isequal(parts.name, ASL_NAME)
          parts.isImage = true;
        end

        % last file in directory
        if (count == numel(files))
          group = [group, parts];
          lastFile = true;
        end
        
        % end of group
        if ((~isempty(currentID)) && (~isequal(currentID, parts.id))) || (lastFile)
          % patients must have 3 files, e.g. t1, pre, post images
          if numel(group) ~= 3
            prev = group(numel(group) - 1);
            error(['patient ', prev.id, ' doesn''t have 3 files.']);
          end
          % we got all the data so process it
          exitCode = process_patient(group);
          if exitCode ~= 0
            return
          end
          group = [];
        end
        
        % add to group
        group = [group, parts];
        currentID = parts.id;
      end
    end
    
    exitCode = 0;
  end

  function exitCode = process_patient(patient)
    t1 = nan;
    preImage = nan;
    postImage = nan;
    id = patient(1).id;

    % pop off t1
    for i = 1:numel(patient)
      if patient(i).isT1
        t1 = patient(i).path;
        patient(i) = [];
        break
      end
    end

    % find the pre/post images
    a = patient(1);
    b = patient(2);
    if a.series > b.series 
      preImage = b.path;
      postImage = a.path;
    else
      preImage = a.path;
      postImage = b.path;
    end

    % output is saved in the following locations:
    %   /outputFolder/patient-id/struc.anat
    %   /outputFolder/patient-id/pre
    %   /outputFolder/patient-id/post

    output = [outputFolder '/' id];
    if ~exist(output,'dir')
        mkdir(output);
    end
    
    % TODO: if we get an error stop the entire program
    exitCode = nii_basil(preImage, t1, '', '', output, false, false, [output '/pre']);
    if exitCode ~= 0
      return
    end
    exitCode = nii_basil(postImage, t1, '', '', output, false, false, [output '/post']);
  end

end