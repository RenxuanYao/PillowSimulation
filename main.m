function pillowSimulation()
    % Clear workspace and close figures
    close all; clc;
    
    % Parameters
    mass = 5;           % Head mass [kg]
    g = 9.81;           % Gravitational acceleration [m/s²]
    force = mass * g;   % Applied force = m*g (Newton's 2nd Law)
    gridSize = [50, 50];% Grid resolution
    stiffness = 5000;   % Material stiffness [N/m²]
    
    % List of available shapes (8 options)
    shapeList = {'flat', 'parabolic', 'wave', 'gaussian', ...
                 'cylindrical', 'conical', 'asymmetric', 'saddle'};
    selectedShape = 'cylindrical'; % Choose shape here
    
    % Create pillow shape
    pillowShape = createPillowShape(gridSize, selectedShape);
    
    % Simulate compression and calculate metrics
    [~, contactArea, pressure] = simulateCompression(pillowShape, stiffness, force);
    avgPressure = force / contactArea;       % Average pressure = Force/Area (Pascal)
    maxPressure = max(pressure(:));          % Peak pressure
    stdPressure = std(pressure(pressure > 0));
    protectionScore = 1 / (maxPressure + stdPressure); % Custom comfort metric
    
    fprintf('\n=== Final Results ===\n');
        fprintf('Contact Area: %.4f m²\n', contactArea);
        fprintf('Average Pressure: %.2f Pa\n', avgPressure);
        fprintf('Max Pressure: %.2f Pa\n', maxPressure);
        fprintf('Pressure Deviation: %.2f Pa\n', stdPressure);
        fprintf('Neck Protection Score: %.6f (higher=better)\n', protectionScore);
    % Visualization
    visualizeResults(pillowShape, pressure, selectedShape);
end

function pillow = createPillowShape(gridSize, type)
    % Creates pillow shapes using geometric formulas.
    % Formula implementations:
    [X, Y] = meshgrid(linspace(-1, 1, gridSize(2)), linspace(-1, 1, gridSize(1)));
    
    switch type
        %-----------------------------------------------
        case 'flat'  % Flat surface
            % Formula: z(x,y) = 0
            pillow = zeros(gridSize);
        
        %-----------------------------------------------
        case 'parabolic'  % Paraboloid shape
            % Formula: z(x,y) = A*(1 - x² - y²)
            A = 0.1; % Amplitude [m]
            pillow = A * (1 - X.^2 - Y.^2);
        
        %-----------------------------------------------
        case 'wave'  % Sinusoidal wave pattern
            % Formula: z(x,y) = A*sin(ω_x*x + ω_y*y)
            A = 0.05; % Amplitude [m]
            Omega_x = 4*pi; Omega_y = 4*pi; % Angular frequencies
            pillow = A * sin(Omega_x*X + Omega_y*Y);
        
        %-----------------------------------------------
        case 'gaussian'  % Gaussian bump
            % Formula: z(x,y) = A*exp(-(x² + y²)/(2σ²))
            A = 0.15; % Amplitude [m]
            Sigma = 0.4;  % Standard deviation
            pillow = A * exp(-(X.^2 + Y.^2)/(2*Sigma^2));
        
        %-----------------------------------------------
        case 'cylindrical'  % Circular ridge
            % Formula: z(x,y) = A*(sqrt(x² + y²) ≤ R)
            A = 0.1; % Height [m]
            R = 0.6; % Radius
            pillow = A * (sqrt(X.^2 + Y.^2) <= R);
        
        %-----------------------------------------------
        case 'conical'  % Linear cone
            % Formula: z(x,y) = A*(1 - sqrt(x² + y²))
            A = 0.12; % Amplitude [m]
            pillow = A * (1 - sqrt(X.^2 + Y.^2));
        
        %-----------------------------------------------
        case 'asymmetric'  % Custom asymmetric shape
            % Formula: z(x,y) = A1*x² + A2*y²
            A1 = 0.08; A2 = 0.12; % Directional amplitudes
            pillow = A1*X.^2 + A2*Y.^2;
        
        %-----------------------------------------------
        case 'saddle'  % Hyperbolic paraboloid
            % Formula: z(x,y) = A*(x² - y²)
            A = 0.1; % Amplitude [m]
            pillow = A * (X.^2 - Y.^2);
    end
end

function [compression, contactArea, pressure] = simulateCompression(pillow, k, F_target)
    % Hooke's Law Simulation: F = k * x (spring force)
    [rows, cols] = size(pillow);
    cellArea = 1 / (rows * cols); % Area per grid cell [m²]
    compression = zeros(size(pillow));
    F_total = 0;
    delta = 0.001; % Compression step [m]
    head_pos = max(pillow(:)) + 0.05; % Initial head position
    
    % Iterative force balance
    while F_total < F_target
        head_pos = head_pos - delta;
        compression = max(0, head_pos - pillow);
        F_total = sum(k * cellArea * compression, 'all'); % Total force
    end
    
    contactArea = sum(compression > 0, 'all') * cellArea;
    pressure = k * compression; % Pressure = Force/Area = k*x (Hooke's Law)
    
    fprintf('Final Head Position: %.4f m\n', head_pos);
    fprintf('Final Force: %.2f N (Target: %.2f N)\n', F_total, F_target);
    fprintf('Final Contact Area: %.4f m²\n\n', contactArea);
end

function visualizeResults(pillow, pressure, shapeName)
    figure;
    surf(pillow, 'EdgeColor', 'none');
    title(['Pillow Shape: ', shapeName]);
    xlabel('X'); ylabel('Y'); zlabel('Height [m]');
    
    figure;
    imagesc(pressure);
    title('Pressure Distribution [Pa]');
    colorbar;
    axis equal tight;
end