import numpy as np

def simulate_shader():
    """
    Simulates the GLSL fragment shader logic in Python to identify the correct image.
    """
    # Shader constants
    shader_size = 20.0
    
    # mat4 points: vec4(r,g,b,a) is used to store 2 vec2 points as (p1.x, p1.y, p2.x, p2.y)
    # The loop goes to i<5, and a mat4 has indices 0-3. Accessing points[4] is
    # out of bounds, which in many GLSL compilers defaults to vec4(0,0,0,0).
    # So we have 10 points in total.
    points_data = [
        [0.1, 0.1], [0.5, 0.5],
        [0.8, 0.8], [0.25, 0.5],
        [0.9, 0.1], [0.5, 0.25],
        [0.1, 0.9], [0.8, 0.2],
        [0.0, 0.0], [0.0, 0.0]  # From the i<5 loop reading points[4]
    ]
    points = np.array(points_data)

    # GLSL functions ported to Python/Numpy
    def dist(a, b):
        """
        Calculates distance, but only if points 'a' and 'b' are in neighboring
        grid cells. 'a' can be an array of UV coordinates.
        """
        A_grid = np.floor(a * shader_size)
        B_grid = np.floor(b * shader_size)
        
        # This is a grid-based optimization. Only calculate distance for nearby points.
        # This logic is equivalent to the GLSL 'if' statement for a whole array 'a'.
        is_neighbor = (np.abs(B_grid[0] - A_grid[:,:,0]) <= 1) & \
                      (np.abs(B_grid[1] - A_grid[:,:,1]) <= 1)
        
        uv = a - b
        euclidean_dist = np.sqrt(uv[:,:,0]**2 + uv[:,:,1]**2)
        
        # If not a neighbor, return a large distance (2.0); otherwise, return the real distance.
        return np.where(is_neighbor, euclidean_dist, 2.0)

    def worley(xy, feature_points):
        """
        Calculates the minimum distance from each pixel in xy to any of the feature_points.
        """
        min_dist = np.full(xy.shape[:2], 2.0)
        for p in feature_points:
            d = dist(xy, p)
            min_dist = np.minimum(min_dist, d)
        return min_dist

    # --- Main fragment() logic ---
    
    # Generate UV coordinates for an output image
    width, height = 60, 30
    u = np.linspace(0, 1, width)
    v = np.linspace(0, 1, height)
    uu, vv = np.meshgrid(u, v)
    uv_grid = np.stack((uu, vv), axis=-1)

    # Calculate W1, W2, W3 as in the shader
    w1 = worley(uv_grid, points)
    # W2 and W3 are calculated in the shader but ultimately discarded
    # w2 = worley(uv_grid + 0.025, points) 
    # w3 = worley(uv_grid - 0.025, points)

    # The final color is vec3(w1), making it grayscale
    final_color_values = w1
    
    # Print the "final equation" for a sample point as requested
    sample_uv = np.array([[[0.5, 0.15]]]) # Using a shape numpy can handle
    sample_w1 = worley(sample_uv, points)[0,0]
    sample_w2 = worley(sample_uv + 0.025, points)[0,0]
    sample_w3 = worley(sample_uv - 0.025, points)[0,0]

    print("--- Shader Calculation for a Sample Pixel ---")
    print(f"For a sample pixel at UV = (0.5, 0.15):")
    print(f"1. Three Worley values are computed:")
    print(f"   W1 = Worley(UV, points) = {sample_w1:.4f}")
    print(f"   W2 = Worley(UV + 0.025, points) = {sample_w2:.4f}")
    print(f"   W3 = Worley(UV - 0.025, points) = {sample_w3:.4f}")
    print(f"2. A temporary vector vec3({sample_w1:.4f}, {sample_w2:.4f}, {sample_w3:.4f}) is created.")
    print(f"3. The '.x' component of this vector is selected: {sample_w1:.4f}")
    print(f"4. The final RGB color is vec3({sample_w1:.4f}), which is ({sample_w1:.4f}, {sample_w1:.4f}, {sample_w1:.4f}).")
    print("This results in a grayscale image based only on the W1 value.")
    print("\n--- ASCII Art representation of the output pattern ---\n")
    
    # Generate and print ASCII art based on the final color values
    chars = ['@', '#', 'S', '%', '?', '*', '+', ';', ':', ',', '.']
    # Normalize the values a bit for better character distribution, clamping high values
    normalized_values = np.clip(final_color_values * 15, 0, len(chars) - 1).astype(int)
    
    for i in range(height):
        print("".join([chars[j] for j in normalized_values[i]]))

simulate_shader()
<<<C>>>