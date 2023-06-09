function varargout = abcd2blocks(varargin)
% Takes a "restricted data" CSV file from the HCP and generates
% a block file that can be used to make permutations in PALM.
%
% Usage:
% [EB,tab] = abcd2blocks(restrfile,blocksfile,dz2sib,ids)
%
% Inputs:
% restrfile  : CSV file from ABCD/NDAR containing at least these fields:
%              subjectid, rel_family_id, rel_relationship, and zygosity fields.
%              Typically, this is the 'acspsw03.txt' file.
% blocksfile : CSV file to be created, with the exchangeability blocks,
%              ready for used with PALM.
% dz2sib     : (Optional) Defines whether dizygotic twins should be
%              treated as ordinary siblings (true), or be a category
%              on its own (false). Default = false.
% ids        : (Optional) A cell array of subject IDs. If supplied, only the
%              subjects with the indicated IDs will be used.
%
% Outputs (if requested):
% EB      : Block definitions, ready for use, in the original order
%           as in the CSV file.
%
% Reference:
% * Winkler AM, Webster MA, Vidaurre D, Nichols TE, Smith SM.
%   Multi-level block permutation. Neuroimage. 2015;123:253-68.
%
% _____________________________________
% Anderson M. Winkler
% NIH/NIMH 
% Dec/2013 (first version, for HCP)
% Apr/2023 (this version, for ABCD)
% http://brainder.org

warning off backtrace
% restrfile = 'acspsw03.txt';
% blocksfile = 'eb.csv';
% dz2sib = false;
% ids = SelIDs;

if nargin >= 1
    restrfile = varargin{1};
else
    error('Input file must be provided, e.g., "acspsw03.txt".');
end
if nargin >= 2
    blocksfile = varargin{2};
else
    blocksfile = [];
end
if nargin >= 3
    dz2sib = varargin{3};
else
    dz2sib = false;
end
if nargin == 4
    ids = varargin{4};
else
    ids = [];
end

% Load the data
acspsw03   = readtable(restrfile);
colnames   = acspsw03.Properties.VariableNames;
acspsw03   = table2cell(acspsw03);

% The rel_relationship column provides the relationship of the participant to his family:
% 0 = single
% 1 = sibling
% 2 = twin
% 3 = triplet

% The genetic_zygosity_status_# provides genetically inferred zygosity status between participant and genetic_paired_subjectid_#
%  1 = monozygotic
%  2 = dizygotic
%  3 = siblings
% -1 = not available

% Locate the relevant column indices from the input file
egid_col    = strcmpi(colnames,'subjectkey');
famid_col   = strcmpi(colnames,'rel_family_id');
rel_gid_col = strcmpi(colnames,'rel_group_id');
zygo1_col   = strcmpi(colnames,'genetic_zygosity_status_1');
zygo2_col   = strcmpi(colnames,'genetic_zygosity_status_2');
zygo3_col   = strcmpi(colnames,'genetic_zygosity_status_3');
zygo4_col   = strcmpi(colnames,'genetic_zygosity_status_4');
event_col   = strcmpi(colnames,'eventname');

% Drop non-baseline rows
idx = ~ strcmpi(acspsw03(:,event_col),'baseline_year_1_arm_1');
acspsw03(idx,:) = [];
Nall = size(acspsw03,1);

% Keep the raw IDs, but use simple numerical IDs
ids_raw   = acspsw03(:,egid_col);
ids_new   = (1:Nall)';

% Merge zygosities:
%   MZ is coded as 1000
%   DZ is coded as 100
%   Singleton or ordinary sib is coded as 10
zygos = cell2mat(acspsw03(:,zygo1_col|zygo2_col|zygo3_col|zygo4_col));
sibtype = 10*ones(size(zygos));
sibtype(zygos == 1) = 1000;
if dz2sib
    sibtype(zygos == 2) = 10;
else
    sibtype(zygos == 2) = 100;
end
sibtype = max(sibtype,[],2);

% Family IDs
famid = cell2mat(acspsw03(:,famid_col));
famid = famid + 1; % there is a family with ID = 0

% Relationships
relgid = cell2mat(acspsw03(:,rel_gid_col));

% Subselect subjects as needed
if numel(ids) > 0
    idx = zeros(size(ids));
    for i = 1:numel(ids)
        idxx = find(strcmpi(ids_raw,ids{i}));
        if not(any(idxx))
            error('Subject %s not found. Remove from the provided set of IDs and re-run.',ids{i});
        end
        idx(i) = idxx;
    end
    ids_new = ids_new(idx);
    famid   = famid(idx,:);
    relgid  = relgid(idx,:);
    sibtype = sibtype(idx,:);
end
N = numel(ids_new);

% If the MZ or DZ twin is missing for a subject, treat them as ordinary sib
twidx = find(sibtype == 1000 | sibtype == 100);
for s = twidx'
    idx = famid == famid(s) & sibtype == sibtype(s);
    if sum(idx) == 1
        sibtype(s) = 10;
    end
end

% Label each family according to their type. The "type" is
% determined by the number and type of siblings.
% Also, get the family size.
famtype = zeros(N,1);
famsize = zeros(N,1);
F = unique(famid);
for f = 1:numel(F)
    fidx = F(f) == famid;
    famtype(fidx) = sum(sibtype(fidx));
    famsize(fidx) = sum(fidx);
end

% For families of size 1, treat MZ or DZ as singletons (missing sib)
sibtype(famsize == 1) = 10;

% Assume rel_group_id is a pregnancy index (a pregnancy can have multiple
% babies). Thus, families that have as many individuals as pregnancies can
% have their individuals freely permuted.
for ft = famid'
    idx = famid == ft;
    tmp = relgid(idx);
    nU = numel(unique(tmp));
    nI = size(tmp,1);
    if nU == nI
        relgid(idx) = 1;
    end
end

% Now the EBs   
B = [-ones(N,1) famtype -famid relgid sibtype ids_new];

% Save as CSV
if ~isempty(blocksfile) && ischar(blocksfile)
    dlmwrite(blocksfile,B,'precision','%d'); %#ok<DLMWT> 
end

% Return EB array to the workspace
if nargout
    varargout{1} = B;
end

% That's it! :-)
