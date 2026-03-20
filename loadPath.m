function pathData = loadPath(folder,blockName,extractMode)
%LOADINFO Load session information from each block directory
%

if ~exist('extractMode', 'var') || isempty(extractMode)
	extractMode = 'struct';
end

pathData = cellfun(@(bl) load(fullfile(folder,'Position',bl, 'pathdata.mat')), blockName, 'un', 0);

if strcmpi(extractMode, 'struct')
	pathData = cell2structure(pathData, blockName);
elseif strcmpi(extractMode, 'cell')
	% pass
end
end

function S = cell2structure(C, fields)

C = cellfun(@(c) c, C, 'un', 0);

args = cell(length(C)*2, 1);
args(1:2:end) = fields;
args(2:2:end) = C;

S = struct(args{:});
end