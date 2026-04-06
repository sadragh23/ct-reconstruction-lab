function A = create_system_matrix(n_detector, angles_deg, n_image, detector_spacing, image_spacing)
% CREATE_SYSTEM_MATRIX Creates system matrix for 2D parallel beam CT
%
%   A = create_system_matrix(n_detector, angles_deg, n_image)
%   A = create_system_matrix(n_detector, angles_deg, n_image, detector_spacing)
%   A = create_system_matrix(n_detector, angles_deg, n_image, detector_spacing, image_spacing)
%
%   Inputs:
%       n_detector        - Number of detector pixels
%       angles_deg        - Array of projection angles in degrees
%       n_image           - Side length of square image (n_image x n_image pixels)
%       detector_spacing  - (optional) Spacing between detector pixels (default: 1.0)
%       image_spacing     - (optional) Spacing between image pixels
%                           (default: auto-sized to fit in reconstruction circle)
%
%   Output:
%       A - System matrix of size (n_angles*n_detector) x (n_image^2)
%           Each row corresponds to one ray (detector pixel at one angle)
%           Each column corresponds to one image pixel
%           A(i,j) = length of ray i passing through image pixel j
%
%   Image setup:
%       - Square image centered at origin
%       - Image domain size depends on n_image and image_spacing
%       - Reconstruction circle diameter = n_detector * detector_spacing
%
%   Geometry:
%       - Parallel beam geometry
%       - Field of view (reconstruction circle) determined by detector width

    % Handle optional arguments
    if nargin < 4 || isempty(detector_spacing)
        detector_spacing = 1.0;  % Default: unit spacing for detector
    end

    if nargin < 5 || isempty(image_spacing)
        % Default: size image so it fits snugly in reconstruction circle
        % Reconstruction circle diameter = detector width
        detector_width = n_detector * detector_spacing;
        % Image diagonal should fit within detector width
        image_spacing = detector_width / (n_image * sqrt(2));
    end

    % Convert angles to radians
    angles_rad = deg2rad(angles_deg);
    n_angles = length(angles_rad);

    % Initialize system matrix
    % Rows: one per ray (n_angles * n_detector rays)
    % Cols: one per image pixel (n_image^2 pixels)
    n_rays = n_angles * n_detector;
    n_pixels = n_image * n_image;
    A = zeros(n_rays, n_pixels);

    % Build system matrix
    ray_idx = 0;
    for angle_idx = 1:n_angles
        theta = angles_rad(angle_idx);

        for det_idx = 1:n_detector
            ray_idx = ray_idx + 1;

            % Get ray origin and direction
            [ray_origin, ray_dir] = get_ray(det_idx, n_detector, theta, ...
                                            detector_spacing);

            % Loop through all image pixels 
            pixel_idx = 0;
            for ix = 1:n_image
                for iy = 1:n_image
                    pixel_idx = pixel_idx + 1;

                    % Get pixel bounds (centered at origin)
                    % ix = column (maps to x), iy = row (maps to y, but inverted)
                    pixel_center_x = (ix - (n_image + 1) / 2) * image_spacing;
                    pixel_center_y = ((n_image + 1) / 2 - iy) * image_spacing;  % Invert y to match Cartesian coords

                    % Pixel is a square with width = image_spacing
                    half_pixel = image_spacing / 2;
                    x_min = pixel_center_x - half_pixel;
                    x_max = pixel_center_x + half_pixel;
                    y_min = pixel_center_y - half_pixel;
                    y_max = pixel_center_y + half_pixel;

                    % Compute intersection length
                    ray_length = ray_square_intersection(ray_origin, ray_dir, ...
                                                         x_min, x_max, y_min, y_max);

                    % Store in system matrix
                    A(ray_idx, pixel_idx) = ray_length;
                end
            end
        end
    end
end


function [ray_origin, ray_dir] = get_ray(det_idx, n_detector, theta, detector_spacing)
% GET_RAY Computes ray line for parallel beam geometry
%
%   For parallel beam at angle theta:
%       - Rays travel in direction: (cos(theta), sin(theta))
%       - Each ray is offset perpendicular to the ray direction
%       - Perpendicular direction: (-sin(theta), cos(theta))
%
%   The ray line passes through a point at perpendicular distance t from origin

    % Perpendicular offset for this detector pixel (centered)
    t = (det_idx - (n_detector + 1) / 2) * detector_spacing;

    % Ray origin: point on ray line closest to origin
    % This point is at perpendicular distance t from origin
    ray_origin = t * [-sin(theta); cos(theta)];

    % Ray direction (parallel beam)
    ray_dir = [cos(theta); sin(theta)];
end


function ray_length = ray_square_intersection(ray_origin, ray_dir, x_min, x_max, y_min, y_max)
% RAY_SQUARE_INTERSECTION Computes length of ray passing through square
%
%   Uses the slab method for ray-box intersection
%   Ray: p(t) = ray_origin + t * ray_dir, for t >= 0
%   Square: [x_min, x_max] x [y_min, y_max]
%
%   Returns: ray_length - length of intersection (0 if no intersection)

    % Small epsilon to avoid division by zero
    epsilon = 1e-10;

    % Initialize intersection parameter range
    t_min = -inf;
    t_max = inf;

    % X-axis slab
    if abs(ray_dir(1)) > epsilon
        tx1 = (x_min - ray_origin(1)) / ray_dir(1);
        tx2 = (x_max - ray_origin(1)) / ray_dir(1);
        t_min = max(t_min, min(tx1, tx2));
        t_max = min(t_max, max(tx1, tx2));
    else
        % Ray parallel to X-axis
        if ray_origin(1) < x_min || ray_origin(1) > x_max
            ray_length = 0;
            return;
        end
    end

    % Y-axis slab
    if abs(ray_dir(2)) > epsilon
        ty1 = (y_min - ray_origin(2)) / ray_dir(2);
        ty2 = (y_max - ray_origin(2)) / ray_dir(2);
        t_min = max(t_min, min(ty1, ty2));
        t_max = min(t_max, max(ty1, ty2));
    else
        % Ray parallel to Y-axis
        if ray_origin(2) < y_min || ray_origin(2) > y_max
            ray_length = 0;
            return;
        end
    end

    % Check if there is a valid intersection
    if t_max < t_min
        ray_length = 0;
        return;
    end

    % Compute intersection length (use full ray line, not just forward direction)
    ray_length = t_max - t_min;
end