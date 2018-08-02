classdef MainClass
    % Main class containing all functions required to perform NR-ICP
        % A non-rigid ICP framework is implemented following 
        % Amberg, B., Romdhani, S. and Vetter, T., 2007, June. 
        % Optimal step nonrigid ICP algorithms for surface registration. 
        % In Computer Vision and Pattern Recognition, 2007. 
        % CVPR'07. IEEE Conference on (pp. 1-8). IEEE.
    
   properties(Access = public)
       
   end
   
   methods
          
      function A = triangulation2adjacency(~, face, vertex)
        %   Compute the adjacency matrix of a given triangulation. 
        %   INPUT: face matrix (mx3); vertex matrix (nx3)
        %   OUTPUT: adjacency matrix
        %   MatLab Toolbox Graph 
        %   Copyright (c) 2009, Gabriel Peyre
          
          f = double(face);

          A = sparse([f(:,1); f(:,1); f(:,2); f(:,2); f(:,3); f(:,3)], ...
                   [f(:,2); f(:,3); f(:,1); f(:,3); f(:,1); f(:,2)], ...
                   1.0);
          A = double(A>0);

          return; 
          
          nvert = max(max(face));
          nface = size(face,1);
          A = spalloc(nvert,nvert,3*nface);

          for i=1:nface          
              for k=1:3
                  kk = mod(k,3)+1;
                  if nargin<2
                      A(face(i,k),face(i,kk)) = 1;
                  else
                      v = vertex(:,face(i,k))-vertex(:,face(i,kk));
                      A(face(i,k),face(i,kk)) = sqrt( sum(v.^2) );    % euclidean distance
                  end
              end
          end 
          % make sure that all edges are symmetric
          A = max(A,A');
      end

      function Ic = adjacency2incidence(~, A)
        %   Converts an adjacency matrix to an incidence matrix 
        %   INPUT: adjacency matrix
        %   OUTPUT: incidence matrix (sparse) 
        %   MatLab Toolbox Graph 
        %   Copyright (c) 2009, Gabriel Peyre
          
          % compute list of edges
          [i,j,s] = find(sparse(A));
          I = find(i<=j);
          i = i(I);
          j = j(I);
          % number of edges
          n = length(i);
          % number of vertices
          nverts = size(A,1);

         % build sparse matrix
         s = [ones(n,1); -ones(n,1)];
         is = [(1:n)'; (1:n)'];
         js = [i(:); j(:)];
         Ic = sparse(is,js,s,n,nverts);
         Ic = Ic';

         % fix self-linking problem (0)
         a = find(i==j);
         if not(isempty(a))
             for t=a'
             Ic(i(t),t) = 1;
             end
         end
      end
      
      function line = createLine3d(~, varargin)
          
          if length(varargin)==1
              error('Wrong number of arguments in ''createLine'' ');
    
          elseif length(varargin)==2
              % 2 input parameters. They can be :
              % - 2 points, then 2 arrays of 1*2 double.
              v1 = varargin{1};
              v2 = varargin{2};
              if size(v1, 2)==3
                  % first input parameter is first point, and second input is the
                  % second point.
                  line = [v1(:,1) v1(:,2) v1(:,3) v2(:,1)-v1(:,1) v2(:,2)-v1(:,2) v2(:,3)-v1(:,3)];    
              elseif size(v1, 2)==1
                  % first parameter is angle with the vertical
                  % second parameter is horizontal angle with Ox axis
                  theta = varargin{1};
                  phi = varargin{2};

                  sit = sin(theta);
                  p0 = zeros([size(theta, 1) 3]);
        
                  line = [p0 cos(phi)*sit sin(phi)*sit cos(theta)];
              end
    
          elseif length(varargin)==3
              % 3 input parameters :
              % first parameter is a point of the line
              % second parameter is angle with the vertical
              % third parameter is horizontal angle with Ox axis
              p0      = varargin{1};
              theta   = varargin{2}; 
              phi     = varargin{3};
    
              if size(p0, 2)~=3
                  error('first argument should be a 3D point');
              end
    
              sit = sin(theta);
              line = [p0 cos(phi)*sit sin(phi)*sit cos(theta)];
    
          elseif length(varargin)==4
              % 4 input parameters :
              p0 = varargin{1};
              dx = varargin{2};
              dy = varargin{3};
              dz = varargin{4};
    
              if size(p0, 2)~=3
                  error('first argument should be a 3D point');
              end
    
              line = [p0 dx dy dz];
              
          elseif length(varargin)==5
              % 5 input parameters :
              % first parameter is distance of lin to origin
              % second parameter is angle with the vertical
              % third parameter is horizontal angle with Ox axis
              x0      = varargin{1};
              y0      = varargin{1};
              z0      = varargin{1};
              theta   = varargin{2}; 
              phi     = varargin{3};
    
              sit = sin(theta);
              line = [x0 y0 z0 cos(phi)*sit sin(phi)*sit cos(theta)];

          elseif length(varargin)==6
              % 6 input parameters, given as separate arguments for x0, y0, z0 and 
              % dx, dy and dz
              % (not very useful, but anyway...)
              x0 = varargin{1};
              y0 = varargin{2};
              z0 = varargin{3};
              dx = varargin{4};
              dy = varargin{5};
              dz = varargin{6};

              line = [x0 y0 z0 dx dy dz];
          end
      end
     
      function [projections] = projectNormals(obj, sourceVertices, Target, normals)
          
          % Get number of source vertices
          nVerticesSource = size(sourceVertices, 1);

          % Pre-allocate space for projections
          projections = zeros(nVerticesSource, 3);

          % Loop over source vertices projecting onto the target surface
          for i=1:nVerticesSource
              
              % Get vertex and normal
              vertex = sourceVertices(i,:);
              normal = normals(i,:);
    
              % Define line in direction normal that passes through vertex
              line = createLine3d(obj, vertex, normal(1), normal(2), normal(3));
    
              % Compute the intersection of the line and the source surface
              intersection = intersectLineMesh3d(line, Target.vertices, Target.faces); 
    
              % If multiple intersections choose the one closest to the source vertex
              if ~isempty(intersection)
                  [~,I] = min(sqrt(sum((intersection - ...
                  repmat(vertex,size(intersection,1),1)).^2, 2)));
                  projections(i,:) = intersection(I,:);
              else
                  % If no intersections just keep the source vertex position
                  projections(i,:) = vertex;
              end 
          end
      end
      
      function [samples] = sampleVerts(~, Mesh, radius )
          
          samples = [];
          vertsLeft = Mesh.vertices;
          itt = 1;
          while size(vertsLeft, 1) > 0
              nVertsLeft = size(vertsLeft, 1);
    
              % pick a sample from remaining points
              vertN = randsample(nVertsLeft, 1);
              vert = vertsLeft(vertN, :);
    
              % Add it to sample set
              samples(itt,:) = vert;
    
              % Remove nearby vertices
              idx = rangesearch(vertsLeft, vert, radius);
              idRemove = idx{1};
              vertsLeft(idRemove, :) = [];

              itt = itt + 1;
          end
       end
      
      function bound = find_bound(~, pts, poly)
           % Finds indices of vertices located at the border of surface
           % INPUT: vertices (nx3), faces (mx3) of target 
           % OUTPUT: rowvector (px1) with p indices refering to vertices at
           % the border 

          TR = triangulation(poly, pts);
          FF = freeBoundary(TR);
          bound = FF(:,1);
      end

      function [ N ] = createD_N(obj, vertsSource, normalsSource )
        %   Create matrix in accordance with formula [8] Amberg et al.  
        %   INPUT: 
        %   vertsSource >> Vertices (nx3)
        %   normalsSource >> Corresponding normals (nx3) 
        %   OUTPUT:
        %   Matrix N (D) >> see equ. [8], nx4n, sparse matrix    
        count = 1; 

        % Initialise empty sparse matrix D and fill it in for loop 
        N = sparse(size(vertsSource, 1), 4*size(vertsSource, 1));
        % Add a column of ones to create homogeneous coordinates
        tempN = ones(size(vertsSource, 1),1);
        homogCoordN = [normalsSource tempN];

        for i=1:size(normalsSource,1)
            N(i,count)= homogCoordN(i, 1);
            N(i,count+1)= homogCoordN(i, 2);
            N(i,count+2)= homogCoordN(i, 3);
            N(i,count+3)= homogCoordN(i, 4);
            count = count + 4;
        end
      end

      function [G_mean] = gaussian_curvature(obj, V, F)
          
          % Computing the Gaussian Curvature using the angle deficit
          % Use the cosine rule since we are looping through 
          % the faces of a triangular mesh
          
          angles = zeros(size(F));
          angles_K = zeros(size(V, 1), 1);
          areas = zeros(size(V, 1), 1);

          for r = 1:size(F, 1)
              
              a = norm ( V( (F(r, 2)) ) - V( (F(r, 3)) ) );
              b = norm ( V( (F(r, 1)) ) - V( (F(r, 3)) ) );
              c = norm ( V( (F(r, 1)) ) - V( (F(r, 2)) ) );

              angles(r, 1) = real( acos( (b^2 + c^2 - a^2) / (2 * b * c) ) );
              angles(r, 2) = real( acos( (a^2 + c^2 - b^2) / (2 * a * c) ) );
              angles(r, 3) = real( acos( (a^2 + b^2 - c^2) / (2 * a * b) ) );

              angles_K( F(r, 1) ) = angles_K( F(r, 1) ) + angles(r, 1);
              angles_K( F(r, 2) ) = angles_K( F(r, 2) ) + angles(r, 2);
              angles_K( F(r, 3) ) = angles_K( F(r, 3) ) + angles(r, 3);

              areas( F(r, 1) ) = areas( F(r, 1) ) + ( ( 0.5 * b * c * sin(angles(r, 1)) ) / 3 );
              areas( F(r, 2) ) = areas( F(r, 2) ) + ( ( 0.5 * b * c * sin(angles(r, 2)) ) / 3 );
              areas( F(r, 3) ) = areas( F(r, 3) ) + ( ( 0.5 * b * c * sin(angles(r, 3)) ) / 3 );  
          end

          G_mean = mean( (2 * pi - angles_K) );

      end
      
      function [ verticesTransformed, X, ErrorVector ] = nricp(obj, Source, Target, Options )
        % Implementation of NR-ICP following Amberg et al., 2007
        %   Inputs:
        %   Source: structure 'Source = struct' with the following fields: 
        %       Source.vertices: nx3 matrix (template vertices)
        %       Source.faces: mx3 matrix (template faces) 
        %       (Optional)Source.normals: nx3 per vertex normals  
        %   Target : structure 'Target = struct' as defined above
        %   Options : Structure with the following fields
        %   gamm:   double, weights differences in the rotational and skew
        %           part of the deformation against the translational part 
        %           Following Amberg: gamma = 1 
        %   epsilon: double, threshold for changes in deformation X.
        %   alphaSet : decreasing vector of (double) stiffness parameters. 
        %           High stiffness parameters force global transformations whereas 
        %           low values allow for local deformations. This describes the
        %           outer loop for decreasing alpha values. (4.1, Page 3)
        %   useNormals : logical (0,1), specifies that surface normals are to be used
        %           to project the source onto the target surface. If this term is
        %           used then the Source input should contain a normals field.
        %   normalWeighting: This will only be enabled if useNormals is enabled as
        %           well. This is an additional weighting factor to wVec (similar
        %           to ignoreBoundary). It sets the weight for a matching point
        %           pair to 0 if the normals vary too much. 
        %   rigidInit : logical (0,1), specifies that rigid ICP should be performed
        %           first before allowing non-rigid and non-global deformations. If
        %           this is set to 0, then a random X will be used. Testing showed
        %           that the alignment still works, however, the result is
        %           different (potentially worse) and it takes longer to run. 
        %   ignoreBoundary: This is an additional weighting factor for how certain 
        %           we are that a pair is actually a match. The weighting is for
        %           wVec in formula 2 (Page 3). The weighting is set to 1, or 0 
        %           if there is no match. 
        %
        %   Outputs:
        %   vertsTransformed : N X 3 vertices of transformed source mesh,
        %   X : (4N) X 3 stacked matrix of transformations.
        %__________________________________________________________________
        
        if ~isfield(Options, 'gamm')
            % Weighting factor in equation [3]
            Options.gamm = 1;
        end
        
        if Options.normalWeighting == 1
            % Set default threshold of 45 degrees used to reject
            % potentially bad matches based on the differences of their
            % normals 
            angle_thresholdD = 45;
            angle_threshold = deg2rad(angle_thresholdD);
        end

        %__________________________________________________________________
        % Get source and target vertices 
        sourceVertices = Source.vertices;
        targetVertices = Target.vertices;
        % Get the number of Source vertices
        nSourceVertices = size(sourceVertices, 1);
        % For this algorithm, it is crucial that the target has more
        % vertices and faces than the source model.
        if nSourceVertices > size(targetVertices, 1)
            disp('* This might not work. The #Target vertices must be higher than the #Template vertices');
        end

        % Optionally get source and target normals 
        if Options.normalWeighting == 1
            normalsSource = Source.normals;
            normalsTarget = Target.normals;
        end

        %__________________________________________________________________
        % Plot source and target surfaces (initial plot)
        % It will be updated for each iteration step (new deformation)  
    
        fig = figure;
        movegui(fig,'northwest');
        PlotTarget = rmfield(Target, 'normals');
        p = patch(PlotTarget, 'facecolor', 'b', 'EdgeColor',  'none', ...
                  'FaceAlpha', 0.5);
        hold on;

        PlotSource = rmfield(Source, 'normals');
        h = patch(PlotSource, 'facecolor', 'r', 'EdgeColor',  'none', ...
            'FaceAlpha', 0.5);
        material dull; light; grid on; xlabel('x'); ylabel('y'); zlabel('z');
        view([60,30]); axis equal; axis manual;
        legend('Target', 'Source', 'Location', 'best')
        drawnow;
    
        %__________________________________________________________________
        % Set matrix G (equation [3] Amberg) 
        diag_elements = [1 1 1 Options.gamm];
        G = diag(diag_elements);

        % Set incidence matrix M (equation [3] Amberg)  
        % A: sparse adjacency matrix (nSource x nSource) with:
        % 1 (connected vertices) and 0 (not connected vertices))
        A = triangulation2adjacency(obj, Source.faces, Source.vertices);
        % M: sparse incidence matrix (nSource x nEdges) with:
        % 1 (connected) and (-1) (not connected) 
        M = adjacency2incidence(obj, A)';

        K = kron(M, G);

        % Set matrix D (equation [8] Amberg)
        D = createD_N(obj, sourceVertices, sourceVertices);

        % Set weights vector (equation [2] Amberg)
        weights = ones(nSourceVertices,1);

        % Get indices of vertices located at the border of the target 
        if Options.ignoreBoundary == 1
            border = find_bound(obj, targetVertices, Target.faces);
        end
        
        %__________________________________________________________________
        % Perform rigid ICP to find initial correspondences X0 (4nx3)
        % (equation [1] Amberg)
        if Options.rigidInit == 1
            disp('* Performing rigid ICP...');
            if Options.ignoreBoundary == 0
                border = 0;
            end
            % MatLab ICP implementation 
            % Martin Kjer and Jakob Wilm, Technical University of Denmark, 2012
            [R, t] = icp(targetVertices', sourceVertices', 50, 'Verbose', true, ...
                         'EdgeRejection', logical(Options.ignoreBoundary), ...
                         'Boundary', border', 'Matching', 'kDtree');

            X = repmat([R'; t'], nSourceVertices, 1);
            % Transform vertices of source
            verticesTransformed = D*X;

            % Update plot: transformed source vertices  
            set(h, 'Vertices', verticesTransformed);
            drawnow;
        else
            % Initialize X0 arbitrarily  
            X = repmat([eye(3); [0 0 0]], nSourceVertices, 1); 
        end
        
        %__________________________________________________________________

        % Perform non-rigid ICP 
        % get number of element in the set of stiffness parameters  
        nAlpha = numel(Options.alphaSet);

        % Outer loop iterates over stiffness parameters alpha.
        % (4.1 Amberg)
        disp('* Performing non-rigid ICP...');
        % Set index to fill residuals vector
        idx = 1;
        for i = 1:nAlpha

            % Update stiffness
            alpha = Options.alphaSet(i);

            % Choose previous X to be very different from X
            prevX = 10*X;

            % Inner loop to compute deformation for given stiffness.
            % Break when consecutive matrices X are similar.
            while norm(X - prevX) >= Options.epsilon 

                % Transform vertices of source
                verticesTransformed = D*X;

                % Update plot for transformed source vertices   
                set(h, 'Vertices', full(verticesTransformed));
                drawnow;

                % Determine closest points on Target for transformed source points  
                targetIdx = knnsearch(targetVertices, verticesTransformed);
                % U: subset of vertices in Target. 
                U = targetVertices(targetIdx,:);

                % Update weights if matches lie on boundary of target
                % ([2], [3] Amberg)
                if Options.ignoreBoundary == 1
                    targetBorder = ismember(targetIdx, border); 
                    % Set weight for vertices at border to 0  
                    weights = ~targetBorder;
                end
                
                % Update weights if normals of matches differ
                if Options.normalWeighting == 1

                    N = createD_N(obj, sourceVertices, normalsSource);
                    % Transform normals 
                    normalsTransformed = N*X;
                    % Find normals of target that have a correspondence 
                    normalsU = normalsTarget(targetIdx,:);
                    % Find angle between the normals 
                    normalDotProd = dot(normalsU, normalsTransformed, 2);
                    normalNormSource = sqrt(sum(normalsTransformed.^2,2));
                    normalNormTarget = sqrt(sum(normalsU.^2,2));
                    angle = acos((normalDotProd)./(normalNormSource .* normalNormTarget));
                    % Check the angle against the angle threshold {0;1}
                    % Update the weightings 
                    weights = weights .* (angle < angle_threshold);
                end

                % Set weight matrix (#source vertices x #source vertices)
                % (equation [7] Amberg)
                W = weights .* speye(nSourceVertices);

                % Set A and B (equation [12] Amberg)
                % Ignore landmarks term  
                A = [alpha * K; W * D];
                B = [zeros(size(M,1)*size(G,1), 3); W * U];

                % Get optimal transformation X
                prevX = X;
                % Solution for (equation [12] Amberg) 
                % Solves the system of linear equation: X = inv(A'A)*A'B
                X = (A' * A) \ (A' * B);
                % Calculate residuals 
                error = norm((A*X -B),'fro')^2;
                ErrorVector(idx) = error;
                idx = idx+1; 
            end
            % End of WHILE (inner loop)
        end
        % End of FOR loop (outer)  

        % Final transformation parameters in X; transform vertices 
        verticesTransformed = D*X;

        % Project along surface normals to Target 
        if Options.useNormals == 1
            disp('* Projecting transformed points onto target along surface normals...');

            % Get template surface normals
            normalsTemplate = Source.normals;

            N = createD_N(obj, sourceVertices, normalsTemplate);
            % Transform the normals of template 
            normalsTransformed = N*X;
            % calculate angle between normals 
            normalsU = normalsTarget(targetIdx,:);
            normalDotProd = dot(normalsU, normalsTransformed, 2);
            normalNormSource = sqrt(sum(normalsTransformed.^2,2));
            normalNormTarget = sqrt(sum(normalsU.^2,2));
            angle = acos((normalDotProd)./(normalNormSource .* normalNormTarget));
            angle_mean = mean(angle(:));
            angle_mean = rad2deg(angle_mean);

            % Project normals to target surface
            verticesTransformed = projectNormals(obj, verticesTransformed, Target, normalsTransformed);
        else
            % Snap template points to closest point on target
            targetIdx = knnsearch(targetVertices, verticesTransformed);
            corTargets = targetVertices(targetIdx,:);
            if Options.ignoreBoundary == 1
                targetBorder = ismember(targetIdx, border);
                weights = ~targetBorder;
            end
            verticesTransformed(weights,:) = corTargets(weights,:);
        end

        % Update plot
        set(h, 'Vertices', verticesTransformed);
        drawnow;
        pause(1);
        % Define global variables for GUI 
        global residuals;
        residuals = ErrorVector;
        
        global normal_mean;
        if Options.useNormals == 1
            normal_mean = angle_mean;
        end
        
        % Compute Gaussian Curvature average difference
        transformed_curvature = gaussian_curvature(obj, verticesTransformed, Source.faces);
        target_curvature = gaussian_curvature(obj, Target.vertices, Target.faces);
        
        global curv_diff;
        curv_diff = abs( transformed_curvature - target_curvature );
        
      end 
         
      function execute(obj)
          
          % Load data
          addpath data;
          
          %Options.gamm = 1;
          
          global src; 
          checkSrc = evalin( 'base', 'exist(''src'',''var'') == 1' );
          if(checkSrc == 1)
              src = strcat('data/', src, '.mat'); load (src);
          else
              error('PLEASE ENSURE TO CHOOSE A SOURCE MESH');
          end
          
          global trgt; 
          checkTrgt = evalin( 'base', 'exist(''trgt'',''var'') == 1' );
          if(checkTrgt == 1)
              trgt = strcat('data/', trgt, '.mat'); load (trgt);
          else
              error('PLEASE ENSURE TO CHOOSE A TARGET MESH');
          end
          
          % Specify that surface normals are available and can be used.
          global useNormals ignoreBoundary rigidInit MaxAlpha MinAlpha StepsAlpha Epsilon;
          
          checkNormals = evalin( 'base', 'exist(''useNormals'',''var'') == 1' );
          if(checkNormals == 1), Options.useNormals = useNormals; Options.normalWeighting = 1;
          else, Options.useNormals = 0; Options.normalWeighting = 0;
          end
          
          checkIgnoreBoundary = evalin( 'base', 'exist(''ignoreBoundary'',''var'') == 1' );
          if(checkIgnoreBoundary == 1), Options.ignoreBoundary = ignoreBoundary;
          else, Options.ignoreBoundary = 0;
          end
          
          checkRigidInit = evalin( 'base', 'exist(''rigidInit'',''var'') == 1' );
          if(checkRigidInit == 1), Options.rigidInit = rigidInit;
          else, Options.rigidInit = 0;
          end
 
          checkMaxAlpha = evalin( 'base', 'exist(''MaxAlpha'',''var'') == 1' );
          checkMinAlpha = evalin( 'base', 'exist(''MinAlpha'',''var'') == 1' );
          checkStepsAlpha = evalin( 'base', 'exist(''StepsAlpha'',''var'') == 1' ); 
          if( (checkMaxAlpha == 0) || (checkMinAlpha == 0) || (checkStepsAlpha == 0))
              disp('Alpha Parameters not completely set, going with default settings ...');
              Options.alphaSet = linspace(100, 10, 20);
          end
          if( ( (checkMaxAlpha == 1) && (checkMinAlpha == 1) && (checkStepsAlpha == 1)) )
              
              if( (isnan(MaxAlpha)) || (isnan(MinAlpha)) || (isnan(StepsAlpha)) )
                  error('PLEASE ENSURE ALPHA PARAMETERS ARE NUMERIC!!!');
              end
              if( (isnumeric(MaxAlpha)) && (isnumeric(MinAlpha)) && (isnumeric(StepsAlpha)) )
                  Options.alphaSet = linspace(MaxAlpha, MinAlpha, StepsAlpha);
              end
          end
          
          checkEpsilon = evalin( 'base', 'exist(''Epsilon'',''var'') == 1' );
          if(checkEpsilon == 1)
              
              if(isnan(Epsilon))
                  error('PLEASE ENSURE EPSILON IS NUMERIC!!!');
              else 
                  Options.epsilon = Epsilon;
              end
          else
              Options.epsilon = 1e-04;
          end
          
          % Perform non-rigid ICP        
          if(trgt == "data/faceTargetMissing.mat")
              nricp(obj, Source, TargetMissing, Options);
          else
              nricp(obj, Source, Target, Options);
          end
      end
           
   end   
end