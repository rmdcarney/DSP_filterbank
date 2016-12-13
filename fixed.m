%
%  fixed:
%     Class implementing fixed point arithmetics.
%     Construct a fixed point variable for a scalar, vector or matrix
%     using f = fixed(bits,value)
% Output:
%       f       - An object representing the fixed point scalar/vector/matrix
% Input:
%       bits    - number of bits including sign bit
%       value   - value
%
%    This fixed point variable can then be used as usual using +, -, *,
%    subscripting and transpose. The result of binary operations uses
%    the lowest resolution of the two operands. Note that division is
%    not implemented, since most fixed point processors don't have it.
%    To convert the fixed point value back to a normal Matlab value,
%    use double(f). The resolution can be retrieved using resolution(f).

%     Mats Bengtsson, 23/6 1999
%     Rickard Stridh, 1/10 1999 value < 0
%     Mats Bengtsson,  9/11 2012. Reimplemented using classdef, 
%     changed implementation so the number of bits include the sign bit
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef fixed
  properties
    bits    % number of bits 
    value   % double value 
  end
  methods
    function f = fixed(bits,value)
      
    % f = fixed(bits,value)
    % Constructor 
    % Output:
    %       f       - An object representing a fixed point number
    % Input:
    %       bits    - number of bits including sign bit
    %       value   - double value 
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      if nargin < 2
        value=0;
        if nargin < 1
          bits=8;
        end
      end
      
      if isa(bits,'fixed')
        f=properties(bits);
          warning('In copy constructor!')
      else % The factor 1-eps makes sure that 1-pow2(-bits) is rounded down.
        f.value=round(double(value) * (1-eps) * pow2(bits-1)) / pow2(bits-1);
        f.bits=bits;
      end
      
      %Overflow check
      if (any(f.value(:) >= 1)||any(f.value(:) < -1))
	maxval=1-pow2(1-bits);
        warning(sprintf('Overflow! Values truncated to %g / -1 !',maxval))
        f.value(f.value>=1)=maxval;
        f.value(f.value< -1)=-1;
      end
    end   
    
    function d = double(f)
    % Convert back to ordinary double
      d=f.value;
    end
    
    function display(f)
    %  Displays an instance of class fixed. Displays only f.value    
      disp(f.value)
    end

    function l = length(f)
    % Length of the vector
      l = length(f.value);
    end
    
    function diff = minus(term1, term2)
      if ~isa(term1,'fixed')
        bits=term2.bits;
      elseif ~isa(term2,'fixed')
        bits=term1.bits;
      else
        bits=min(term1.bits,term2.bits);
      end
      diff = fixed(bits,double(term1)-double(term2));
    end
  
    function prod = mtimes(term1, term2)
      if ~isa(term1,'fixed')
        bits=term2.bits;
      elseif ~isa(term2,'fixed')
        bits=term1.bits;
      else
        bits=min(term1.bits,term2.bits);
      end
      prod = fixed(bits,double(term1)*double(term2));
    end
  
    function prod = times(term1, term2)
      if ~isa(term1,'fixed')
        bits=term2.bits;
      elseif ~isa(term2,'fixed')
        bits=term1.bits;
      else
        bits=min(term1.bits,term2.bits);
      end
      prod = fixed(bits,double(term1).*double(term2));
    end
    
    function sum = plus(term1, term2)
      if ~isa(term1,'fixed')
        bits=term2.bits;
      elseif ~isa(term2,'fixed')
        bits=term1.bits;
      else
        bits=min(term1.bits,term2.bits);
      end
      sum = fixed(bits,double(term1)+double(term2));
    end

    function s = size(f,varargin)
    % Size of the vector/matrix
      if nargin > 1
        s = size(f.value,varargin{2:end});
      else
        s = size(f.value);
      end
    end
    
    function A = subsasgn(A, S, B)

    % A = subsasgn(A, S, B)
    % Subscript Designator for class fixed
    % Assigns an element in a class fixed instance
    % A = SUBSASGN(A,S,B) is called for the syntax A(I)=B, A{I}=B, or
    % A.I=B when A is an object.  S is a structure array with the fields:
    %     type -- string containing '()', '{}', or '.' specifying the
    %             subscript type.
    %     subs -- Cell array or string containing the actual subscripts. 

      C = fixed(A.bits,B);
      A.value=subsasgn(A.value,S,C.value); 
    end
    
    function B = subsref(A, S)

    % B = subsref(A, S)
    % Subscripted Reference
    % B = SUBSREF(A,S) is called for the syntax A(I), A{I}, or A.I
    % when A is an object.  S is a structure array with the fields:
    %     type -- string containing '()', '{}', or '.' specifying the
    %             subscript type.
    %     subs -- Cell array or string containing the actual subscripts.

      B = fixed(A.bits, subsref(A.value,S));
    end
    
    function t=transpose(f)
    %  t=f.'
      t=fixed(f.bits,f.value.');
    end
  
    function t=ctranspose(f)
    %  t=f'
      t=fixed(f.bits,f.value');
    end
    
    function b=resolution(f)
    % Return the number of bits resolution.
        b=f.bits;
    end
    
    function m = uminus(f)
      m = fixed(f.bits,-f.value);
    end
  end
end