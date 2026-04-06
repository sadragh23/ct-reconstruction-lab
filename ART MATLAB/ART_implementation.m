
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  create a test phantom
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear
clc
% image resolution
N = 100; 

% Create the phantom matrix O (all zeros)
O = zeros(N, N);

% Add a recognizable pattern (two squares)
% This assigns the value 1.0 to a specific range of rows and columns
O(15:35, 15:35) = 1.0; 

O(40:45, 10:20) = 0.5;


% Store the vector-reformatted version as VO
VO = O(:); 

% Visualize the phantom (true µ)
imshow(O);
title('Original Phantom (Matrix O)');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting up the 'Scanner'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


n_detector = round(N * sqrt(2));       % How many 'sensors' the scanner has
                       % (should be enogh to cover the 100 x 100 phantom even diagonal)
                       % diagonal is suqare root of 2
                    
angles_deg = 0:2:178;  % The angles the scanner rotates (0, 2, 4... 178)

% Now we CALL the function to build the matrix D
% D is the system matrix: each row is one X-ray path, each column is one image pixel.
% This uses the variables we defined (n_detector, angles_deg, N)
D = create_system_matrix(n_detector, angles_deg, N);

%%%%%%%%%%%%%%%%%%
% The Forward Projection (P = D*VO) ---
% This simulates shooting the X-rays through our phantom VO
P = D * VO; % Forward Projection vector

% Visualize the result (The Sinogram)
sinogram = reshape(P, n_detector, length(angles_deg));
figure; 
imagesc(angles_deg, 1:n_detector, sinogram);
colormap(gray); 
colorbar;
title('The Sinogram (Our "Measured" Data)');
xlabel('Angle (Degrees)');
ylabel('Detector Pixel Index');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%        ART
%%%%%%%%%%%%%%%%%%%%%%%%%%


% Prepare the Blank Canvas

% Create an empty vector (all zeros) for the reconstruction
VART = zeros(size(VO)); 

% Define how many times we want to loop through the whole dataset

num_iterations = 20; 

% Get the total number of rays (number of rows in matrix D)
n_rays = size(D, 1);

fprintf('Ready to start ART with %d rays and %d iterations.\n', n_rays, num_iterations);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% The ART Reconstruction

for iter = 1:num_iterations
    for i = 1:n_rays
        %  Get the current ray (row i of matrix D)
        ai = D(i, :); 
        
        %  Calculate the squared norm (sum of squares of distances in the row)
        % This is the denominator ||ai||^2 in the formula
        norm_sq = sum(ai.^2);
        
        %  Only update if the ray actually passes through the image
        if norm_sq > 0
            % Calculate the update: (Real Measurement - Predicted) / Norm
            % Then multiply by the ray weights (ai)
            
            update = ((P(i) - ai * VART) / norm_sq) * ai';
            
            % Update our reconstruction vector
            VART = VART + update;
        end
    end
    
    % Print progress
    fprintf('Completed iteration %d of %d\n', iter, num_iterations);
end

% Display the Result
% Reshape the vector back into an N x N image
reconstruction = reshape(VART, N, N);

figure;
subplot(1,2,1); imshow(O); title('Original');
subplot(1,2,2); imshow(reconstruction, []); title(['ART Reconstruction', newline, num2str(iter), ' iterations']);

% error between image and reconstructed image

VSUB = VO - VART;
reconstruction_error = reshape(VSUB, N, N);

figure;
imshow(reconstruction_error); title(['Reconstruction Error', newline, 'after ' num2str(iter), ' iterations']);
