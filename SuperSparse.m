% http://en.wikipedia.org/wiki/Automatically_Tuned_Linear_Algebra_Software
% http://en.wikipedia.org/wiki/Basic_Linear_Algebra_Subprograms
% "Highly optimized implementations have been developed by hardware vendors such as Intel and AMD"



% SUPERSPARSE               sparse arrays of any size and any type
%
% SUPERSPARSE extends Matlab's SPARSE by working aroung some of its
% limitations.
%
% Data held by a SUPERSPARSE can:
%   - be 1, 2, or more than 2 dimensional
%   - contain any data type supported in Matlab (those listed in
%     'help class', as well as any Java classes or user-defined types
%   - exceed the limit on the maximum addressable element implied by
%     32-bit or 64-bit Matlab.
%
% The only limitation in SUPERSPARSE is set by the amount of available
% memory.
%
% SUPERSPARSE(S) will initialize a SUPERSPARSE array of type 'double' 
% with dimensions defined in the vector S. All elements will initially 
% be zero.
%
% SUPERSPARSE(S, 'class') will initialize a SUPERSPARSE array of type 
% 'class'.The type 'class' must be one of those listed in 'help class', 
% or a valid Java class or user-defined class on the Matlab path. All 
% elements will initially be zero (when possible) or empty (when 
% appropriate). Note that like the 'double' sparse, the zeros or 
% empties will not be stored explicitly, but only generated upon 
% referencing empty elements.
%
% SUPERSPARSE(S,I,D) initializes a SUPERSPARSE with dimensions defined 
% in the vector S. The data contained in D will be assigned to indices
% contained in matrix/vector I. The indices I may be linear indices 
% (1-D vector for a non-1-D size S) or subscripts (N-D matrix for an 
% N-D size S). 
%
% S and/or I may be empty ([]), in which case the maximum indices in 
% each dimension of I will determine the size of the SUPERSPARSE. Also 
% I may be left empty, in which case the SUPERSPARSE will have the 
% same size as D.
%
% The type of the SUPERSPARSE array will be identical to that of D, 
% even if D is a cell array containing multiple types. That is, a
% SUPERSPARSE may contain multiple classes. In this case, the implicit 
% 'zeros' will be equal to the empty array ([]). 
%
% SUPERSPARSE(S,I,D, nzMax) will allocate memory for a total of nzMax 
% non-empty elements. Pass an empty array ([]) for I and D when it is 
% known how many elements will be non-empty, but no initial data needs 
% to be set. Note that unlike SPARSE, it is not possible in SUPERSPARSE 
% to allocate more non-empty elements after initializing it with nzMax. 
% To create a dynamic SUPERSPARSE that does not have this restriction, 
% simply leave out the nzMax argument.
%
% SUPERSPARSE(S,I,D,nzMax, zero) will use 'zero' for the implicitly 
% stored zeros. This is useful for example when creating a SUPERSPARSE 
% of a user-defined class, which has a non-trivial way to initialize 
% its equivalent of zero. Leave preceeding arguments empty ([]) when 
% only the last argmument is needed. 
%
%
% EXAMPLE:  (sparse of type 'uint8' of size 2-by-2-by-3)
%
%
% EXAMPLE:  (exceed 64-bit boundary)
%
%
% EXAMPLE:  (sparse of char with custom zero)
%  
%    >> C = SuperSparse([3 3], [1 2],[1 1], 'PQ', [], 'ø')
%
%
%
% EXAMPLE:  (sparse of various types)
%
%
% See also sparse, class, cell.

% Author
%{
    Rody Oldenhuis
    oldenhuis@gmail.com
    LuxSpace Sarl.
%}

% Changelog
%{
    2013/November/06 (Rody Oldenhuis)
        Initial version.
%}
classdef SuperSparse < handle
    
    properties (Access = private, Hidden)
        
        data   % The actual data. Its dimensions and shape, class, etc.
        % depend on the value of the 'type' property.
        
        dims   % Dimensions of this SuperSparse.
        
        inds   % LINEAR indices of non-zero elements
        
        type   % indicator to remember which type of SPARSE this is.
        
        nzMax  % What is the maximum number of non-empty elements we're going
        % to have to store?
        
        maxSz  % maximum number of elements allowed in a matrix
        % on this version of MATLAB.
        
        zero   % Store one 'zero' for the particular class of the SUPERSPARSE.
        % This makes it easier to generate 'zeros' when requested.
    end
    
    methods
        
        %% Constructor
        
        function this = SuperSparse(varargin)
            
            % get max. variable size on this platform
            [DUMMY_, this.maxSz] = computer; %#ok<ASGLU>
            clear DUMMY_ % NOTE: for backwards compatibility
            
            % Parse input arguments
            switch nargin                
                case 0
                    % empty sparse; do nothing.
                    
                case 1
                    % Check for copy constructor
                    if ~isa(varargin{1}, 'SuperSparse')                        
                        this.dims = checkDims(varargin{1});                        
                        this.type = 'double';
                        this.zero = 0;
                        
                    % This is a copy constructor
                    else
                        this.data  = varargin{1}.data;
                        this.dims  = varargin{1}.dims;
                        this.inds  = varargin{1}.inds;
                        this.type  = varargin{1}.type;
                        this.nzMax = varargin{1}.nzMax;
                        this.maxSz = varargin{1}.maxSz;
                        this.zero  = varargin{1}.zero;
                    end
                    
                case 2
                    this.dims = checkDims(varargin{1});
                    this.type = checkType(varargin{2});                    
                    %this.zero =
                    
                case 3
                    this.dims = checkDims(varargin{1});
                    this.inds = checkInds(varargin{2});
                    this.data = checkData(varargin{3});                    
                    this.type = class(this.data);
                    %this.zero =
                    
                case 4
                    this.dims  = checkDims (varargin{1});
                    this.inds  = checkInds (varargin{2});
                    this.data  = checkData (varargin{3});
                    this.nzMax = checkNzMax(varargin{4});
                    this.type  = class(this.data);
                    %this.zero =
                    
                case 5
                    this.dims  = checkDims (varargin{1});
                    this.inds  = checkInds (varargin{2});
                    this.data  = checkData (varargin{3});
                    this.nzMax = checkNzMax(varargin{4});
                    this.type  = class(this.data);
                    this.zero  = varargin{5};
                    
                otherwise
                    error('SuperSparse:SuperSparse:incorrect_argument_count',...
                        'Invalid number of arguments. Type ''help SuperSparse'' for more information.');
            end
            
            % Input checkers
            
            function dims = checkDims(dims)
                if isnumeric(dims)
                    if ~isempty(dims) && (~isvector(dims) || ~all(isfinite(dims)))
                        error('SuperSparse:SuperSparse:checkDims:invalid_dimensions',...
                            'Invalid dimensions specified: ''%s''.', num2str(dims));
                    end
                else
                    error('SuperSparse:SuperSparse:checkDims:invalid_class',...
                        'Function is not defined for dimensions of class ''%s''.', class(dims));
                end
            end
            
            function type = checkType(type)
                if ~ischar(type)
                    error('SuperSparse:SuperSparse:checkType:invalid_type',...
                        'The type of the SuperSparse must be specified with a string.');
                end
                if ~exist(type, 'class')
                    error('SuperSparse:SuperSparse:checkType:invalid_class',...
                        'The specified class ''%s'' cannot be found on the MATLAB path.', type);
                end
            end
            
            function inds = checkInds(inds)
                % TODO
            end
            
            function data = checkData(data)
                % TODO
            end
            
            function nzMax = checkNzMax(nzMax)
                % TODO
            end
            
            % Check consistency of SuperSparse
            % TODO: check datatypes, sizes, nzMax exceeds maxSz, etc.
            
        end
        
        %% Basic functionality
        
        
        % SUBSREF, SUBSASGN and END
        % ----------------------------
        
        
        function out = end(this,k,n) %#ok<INUSD>
            out = this.dims(k); end              
        
        function newSparse = subsref(this, S)
            
            % Basic checks
            switch S.type
                case '{}'
                    error('SuperSparse:subsref:invalid_cell_reference',...
                        'Cell contents reference from a non-cell array object.');
                case  '.'
                    error('SuperSparse:subsref:invalid_structure_reference',...
                        'Attempt to reference field of non-structure array.');
            end    
            
            colons = cellfun('isclass', S.subs, 'char');
%             if any(colons)
%                 S.subs(colons) = 
%             end
            
            if ~all(cellfun(@isnumeric, S.subs))
                ME = MException('SuperSparse:subsref:bad_subscript',...
                    'Indices to a SuperSparse must be of numeric class.');
                throwAsCaller(ME);
            end   
            
            if any(cellfun(@(x)any(x<=0), S.subs))
                ME = MException('SuperSparse:subsref:bad_subscript',...
                    'Subscript indices must either be real positive integers or logicals.');
                throwAsCaller(ME);
            end
            
            % Do we have empty ranges?
            empties = cellfun('isempty', S.subs);                
            if find(~empties,1,'last') > numel(this.dims) || ...
                    any(cellfun(@(x,y) any(x>y), S.subs(~empties), num2cell(this.dims(~empties))))
                ME = MException('SuperSparse:subsref:out_of_bounds',...
                    'Index exceeds SuperSparse dimensions.');
                throwAsCaller(ME);
            end              
            if any(empties)  
                newDims = this.dims;                    
                newDims( empties) = 0; 
                newDims(~empties) = cellfun(@numel, S.subs(~empties));                
                newSparse = SuperSparse(newDims, this.type);                
                return;
            end
            
            % Normal indexing
            [I{1:numel(S.subs)}] = ndgrid(S.subs{:});            
            linInds = sub2ind(this.dims, I{:});
            
            mapInds = ismember(linInds, this.inds)
            
            
            
            
        end
        
        function subsasgn(this, S, vals)
            % TODO
        end
        
        
        % Mimick DISP of SPARSE
        function disp(this)  
            
            % SuperSparse is empty
            if isempty(this)                   
                if ischar(this.type) || numel(this.type)==1
                    if isempty(this.dims)
                        sz = '0-by-0';
                    else
                        sz = regexprep(num2str(this.dims), '\s*', '-by-');
                    end
                    disp(['   All zero ''' this.type ''' SuperSparse: ' sz]);
                elseif isempty(this.type)
                    disp('   Empty SuperSparse.');
                end
                
                return;
            end
            
            % SuperSparse is NOT empty
            
        end
        
        % class, isempty, any, all, reshape, size, numel, ndims, length,
        % permute, squeeze, shiftdim, flipdim, flipud, fliplr, circshift
        % ----------------------------
        
        function str = class(this)
            str = this.type; end
        
        function y = isempty(this)
            y = isempty(this.dims) || any(this.dims==0); end
        
        function y = any(this)
            y = numel(this.data)~=0; end
        
        function y = all(this)
            y = numel(this.data)==prod(this.dims); end
        
        function n = numel(this)
            n = prod(this.dims); end
        
        function n = ndims(this)
            n = numel(this.dims); end
        
        function L = length(this)            
            L = ~isempty(this) * max(size(this)); end
        
        function varargout = size(this, varargin)            
            if nargin == 1
                sz = this.dims;
            elseif nargin == 2
                D = varargin{1};
                if D <= numel(this.dims)
                    sz = this.dims(D);
                else
                    sz = 1;
                end
            else
                error('SuperSparse:size:argin_count',...
                    'Too many input arguments.');
            end
            
            if nargout <= 1
                varargout{1} = sz;
            elseif nargout <= numel(sz)
                varargout = num2cell(sz);
            else
                error('SuperSparse:size:argout_count',...
                    'Too many output arguments.');                
            end            
        end
        
        function newSparse = reshape(this, varargin)
            % TODO
        end
        
        function newSparse = permute(this, varargin)
            % TODO
        end
        
        function newSparse = squeeze(this, varargin)
            % TODO
        end
        
        function [newSparse, shifts] = shiftdim(this, varargin)
            % TODO
        end
        
        function newSparse = flipud(this)
        end
        
        function newSparse = fliplr(this)
        end
        
        function newSparse = flipdim(this, dim)
        end
        
        function newSparse = circshift(this, shiftsize)
        end
        
        
        % Type & attribute checkers
        % ----------------------------
            
        function y = isfloat(this) 
            y = isfloat(this.data) && isfloat(this.zero); end
        
        function y = isinteger(this)
            y = isinteger(this.data) && isinteger(this.zero); end
        
        function y = islogical(this)
            y = islogical(this.data) && islogical(this.zero); end
        
        function y = ischar(this)
            y = ischar(this.data) && ischar(this.zero); end
        
        function y = isa(this, className)
            y = isa(this.data, className) && isa(this.zero, className); end
        
        function y = issparse(this) %#ok<MANU>
            y = true; end
        
        function y = isnumeric(this)
            y = isnumeric(this.data) && isnumeric(this.zero); end
        
%         function y = isnan(this)
%         end
%         
%         function y = isinf(this)
%         end
%         
%         function y = isfinite(this)
%         end

    
        % Type casting 
        % ----------------------------
        
        % TODO
        
        
        % Logical operators
        % ----------------------------
        
        function newSparse = eq(this, other)
        end
        
        function newSparse = ne(this, other)
        end
        
        function newSparse = lt(this, other)
        end
        
        function newSparse = gt(this, other)
        end
        
        function newSparse = le(this, other)
        end
        
        function newSparse = ge(this, other)
        end
        
        
        
        % Math operators
        % ----------------------------
        
        % Unary plus and unary minus: delegate call to class of data
        function this = uplus(this)
            this.data = +this.data; end
        
        function this = uminus(this)
            this.data = -this.data; end
        
        
        function newSparse = plus(this, other)
            try
                this.checkOther(other);
            catch ME
                throwAsCaller(ME); 
            end
            
        end
        
        function newSparse = minus(this, other)
        end
        
        function newSparse = times(this, other)
        end
        
        function newSparse = power(this, other)
        end
        
        function newSparse = mpower(this, other)
        end
        
        function newSparse = mtimes(this, other)
        end
        
        
        
        function newSparse = ldivide(this, other)
        end
        
        function newSparse = rdivide(this, other)
        end
        
        function newSparse = mldivide(this, other)
        end
        
        function newSparse = mrdivide(this, other)
        end
        
        function newSparse = idivide(this, other)
        end
        
        % Sum, diff, bsxfun, repmat, min, max, sort, median, mean
        % ----------------------------
        
        function newSparse = sort(this, varargin)
        end
        
        function newSparse = min(this, varargin)
        end
        
        function newSparse = max(this, varargin)
        end
        
        function newSparse = mean(this, varargin)
        end
        
        function newSparse = median(this, varargin)
        end
        
        function newSparse = sum(this, varargin)
        end
        
        function newSparse = diff(this, varargin)
        end
        
        function newSparse = bsxfun(this, varargin)
        end
        
        function newSparse = repmat(this, varargin)
        end
        
        
    end
        
    %% Internal affairs
    
    methods (Access = protected, Hidden)
    end
    
    methods (Access = private, Hidden)
        
        % Check the 'other' argument -- is it suitable for use in 
        % one of the mathematical operators?
        function checkOther(this, other) %#ok<MANU>
           if ~isnumeric(other) 
               stack = dbstack('-completenames');    
               caller = builtin('_brace', regexp(stack(2).name, '\.(\w*)$', 'tokens'),1);
               ME = MException('SuperSparse:checkOther:invalid_type',...
                   ['Undefined method ''' caller{:} ''' for ''SuperSparse'' and ''' class(other) '''.']);
               throwAsCaller(ME);
           end
        end
        
    end
    
end



