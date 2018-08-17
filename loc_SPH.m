classdef loc_SPH < handle
    %LOC_SPH Summary of this class goes here
    %   Detailed explanation goes here
    
    properties%( SetAccess = protected )
        thr         % Threshold radius
        dim         % Dimensionality of the geometry
        phySize     % The physical dimensions of system
        voxSize     % Number of voxels per dimension
        voxPts      % Points within voxel
        ptsCor      % Point actual coordinates
        ptsLoc      % Points location
        ptsVox      % Points in which voxel
        ptsNum      % Points number
        conList     % List of which particle interacts with which
        conDist     % List of connection distances
        conRHat     % Connection hat vectors
    end
    properties%( Access = protected )
        voxNum      % Voxel count; prod(voxSize)
        voxLink     % Which voxel neighbours which (in linear index)
        voxPtsNum   % Number of points within voxel
        voxL2S      % Function to convert linear index to subscript
        voxS2L      % The reverse function
        voxLoc      % Location on the voxels
        disPos      % Positive displacement vectors
        disNum      % Number of positive displacement vectors
    end

    methods
        function o = loc_SPH(rad,geo)
            %LOC_SPH Construct an instance of this class

            % INPUT CHECKING
            % Rotate geo if not horizontal
            if size(geo,1) ~= 1
                geo = geo';
            end
            % Check if radius makes sense
            if ( ~isscalar(rad) ) || ( ~isreal(rad) ) || ( rad(1) <= 0 )
                error('Radius needs to be a >0 real scalar');
            end
            % Check if geometry is divisible;
            if ~all( mod(geo,rad) == 0 )
                error( 'The geometry must be divisible with grid size!' );
            end
            
            % Data entry
            o.thr = rad;
            o.phySize = geo;
            o.dim = length( geo );
            o.voxSize = geo / rad;
            o.voxNum = prod( o.voxSize );

            % Write down converters
            % Do note that the converter function require column input
            o.voxS2L = @(u) 1+sum((u-1).*[1,cumprod(o.voxSize(1:(end-1)))],2);
            % Linear to subscript is a matrix
            o.voxL2S = cell(1,o.dim);
            [o.voxL2S{:}] = ind2sub( o.voxSize , (1:o.voxNum)' );
            o.voxL2S = cell2mat( o.voxL2S );

            % Calculate positive vectors
            o.disNum = ( 3^o.dim - 1 ) / 2;
            % Initiate for N = 1
            o.disPos = 1;
            sto = (-1:1)';
            for d = 2:o.dim
                o.disPos = [ o.disPos, zeros( size( o.disPos, 1 ), 1 ); ...
                    sto, ones( size( sto, 1 ), 1) ];
                sto = [ kron( ones( 3, 1 ), sto ), ...
                    kron( (-1:1)', ones( size( sto, 1 ), 1 ) ) ];
            end

            % Get movement matrix
            o.voxLink = zeros( o.voxNum, o.disNum );
            for l = 1:o.disNum
                o.voxLink(:,l) = o.voxS2L( ...
                    1 + mod( o.voxL2S + o.disPos(l,:) - 1 , o.voxSize ) );
            end
        end
        
        function insertPoints(o,pts)
            %INSERTPOINTS: Insert points into the instance
            
            % Input check
            if ( length(size(pts)) ~= 2 ) || ( ~any( size(pts) == o.dim ) )
                error( 'Invalid points dimension' );
            end
            if size(pts,2) ~= o.dim
                pts = pts';
            end
            if any( max(pts,[],1) >= (o.phySize) ) || any( min(pts,[],1) < 0 )
                warning('Points out of bound, wrapping around . . .');
                pts = mod( pts, o.phySize );
            end

            % Insert data
            o.ptsCor = pts;
            o.ptsLoc = mod( pts, o.thr );
            o.ptsNum = size(pts,1);
            o.ptsVox = o.voxS2L( 1 + floor( pts / o.thr ) );

            % Store voxel data
            o.voxPts = cell( o.voxNum, 1 );
            o.voxLoc = cell( o.voxNum, 1 );
            o.voxPtsNum = zeros( o.voxNum, 1 );
            for v = 1:o.voxNum
                o.voxPts{v} = find( o.ptsVox == v )';
                o.voxLoc{v} = o.ptsLoc( o.voxPts{v}, : );
                o.voxPtsNum(v) = length( o.voxPts{v} );
            end
        end
        
        function movePoints(o,del)
            %MOVEPOINTS: Move the points by the specified displacement vector
            o.ptsLoc = o.ptsLoc + del;

            % FILL IN LATER
        end
        
        function D = pairwiseDist(o)
            %PAIRWISEDIST Return pairwise distances under length threshold
            % Allocate distance matrix
            allocPts = sum( o.voxPtsNum .^ 2 + ...
                sum( o.voxPtsNum .* o.voxPtsNum(o.voxLink) , 2 ) );
            D2 = sparse( [], [], [], o.ptsNum, o.ptsNum, allocPts );
            
            for v = 1:o.voxNum
                % Get the points within this voxel
                pt1 = o.voxPts{v};
                co1 = o.ptsLoc(pt1,:);
                % Do points within the voxel.
                sq1 = dot(co1,co1,2);
                D2(pt1,pt1) = triu( sq1 + sq1' - 2*(co1*(co1')) , 1 );
                % Do crosstalk between voxels
                for l = 1:o.disNum
                    pt2 = o.voxPts{ o.voxLink(v,l) };
                    co2 = o.ptsLoc(pt2,:) + o.thr * o.disPos(l,:);
                    D2(pt1,pt2) = sq1 + dot(co2,co2,2)' - 2*(co1*(co2'));
                end
            end
            
            D2( D2 >= (o.thr.^2) ) = 0;
            D2 = sqrt( D2 );
            D = sqrt( triu( D2, 1 ) + tril( D2, -1 )' );
            
            % Write results to struct
            o.conList = zeros( nnz(D), 2 );
            [o.conList(:,1),o.conList(:,2),o.conDist] = find( D );
            o.conRHat = o.ptsCor(o.conList(:,2),:) - o.ptsCor(o.conList(:,1),:);
            o.conRHat = o.conRHat ./ sum( o.conRHat.^2, 2 );
            
        end
        
        function D = pairwiseDistBrute(o)
            %PAIRWISEDIST Calculate all to all distance with thresholding
            D = squareform(pdist(o.ptsLoc));
            D( D > o.thr ) = 0;
            D = sparse(triu(D,1));
            
            % Write results to struct
            o.conList = zeros( nnz(D), 2 );
            [o.conList(:,1),o.conList(:,2),o.conDist] = find( D );
            o.conRHat = o.ptsCor(o.conList(:,2),:) - o.ptsCor(o.conList(:,1),:);
            o.conRHat = o.conRHat ./ sum( o.conRHat.^2, 2 );
        end
        
        function plot(o)
            %PLOT Plot points and the detected connections
            
            % Check dimensions for 2d
            if o.dim ~= 2
                error('Plotting done only in 2D');
            end
            
            % Calculate distances
            o.pairwiseDistBrute;
            
            % Get points and ghost points
            col = hsv2rgb( [ linspace(0,1,o.ptsNum)', ones(o.ptsNum,2) ] );
            gho_cor = o.ptsCor + o.phySize .* ( o.ptsCor < o.thr ) ...
                - o.phySize .* ( o.ptsCor >= ( o.phySize - o.thr ) );
            gho = any( ( o.ptsCor < 0 ) | ( o.ptsCor > o.phySize ) , 2 );
            gho_col = col( gho, : );
            
            % Get plotting coordinates
            x1 = o.ptsCor( o.conList(:,1), 1 )';
            y1 = o.ptsCor( o.conList(:,2), 2 )';
            x2 = x1 + (o.conDist') .* (o.conRHat(:,1)');
            y2 = y1 + (o.conDist') .* (o.conRHat(:,2)');
            x = [x1;x2];
            y = [y1;y2];
            
            p_x = o.ptsCor(:,1)';
            p_y = o.ptsCor(:,2)';
            g_x = gho_cor(:,1)';
            g_y = gho_cor(:,2)';
            
            % Plot the thing
            clf
            scatter( p_x, p_y, 20, col, 'filled' );
            hold on
            %scatter( g_x, g_y, 20, gho_col );
            
            % Draw lines
            plot(x,y,'Color',[0,0,0,.3]);
            
            % Axis stuff
            axis([-o.thr,o.phySize(1)+o.thr,-o.thr,o.phySize(2)+o.thr]);
        end
    end
end

