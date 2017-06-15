function p=plotGridPosition_new(epos,nchans,ncols)
%PURPOSE: GIVEN ELECTRODE POSITION, GIVES POSITION OF SUBPLOT IN 16 X 16
%GRID. REDUCES GRAY SPACE.
%INPUTS: EPOS= ELECTRODE POSITION (1-256)

if ~exist('nchans','var') || isempty(nchans)
    nchans = 256;
end

if ~exist('ncols','var') || isempty(ncols)
    ncols = 16;
end

gap = .01;

nrows=ceil(nchans/ncols);

row = nrows - ceil(epos/ncols);
col = rem(epos - 1,ncols);


p(2) = gap+row/nrows*(1-2*gap);
p(1) = gap+col/ncols*(1-2*gap);
p(3) = (1-gap*(ncols+1))/ncols;
p(4) = (1-gap*(nrows+1))/nrows;
hold on