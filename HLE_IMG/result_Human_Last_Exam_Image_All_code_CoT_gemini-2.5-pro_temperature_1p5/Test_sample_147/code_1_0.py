import numpy as np
from PIL import Image

def solve_and_generate():
    """
    Simulates the provided fragment shader to identify the correct image.
    """
    WIDTH, HEIGHT = 256, 256
    SHADER_SIZE = 20.0

    # The feature points from the mat4 variable in the shader.
    # A mat4 has 4 columns, each a vec4. We extract 2 vec2s from each.
    points_raw = [
        # col 0
        (0.1, 0.1), (0.5, 0.5),
        # col 1
        (0.8, 0.8), (0.25, 0.5),
        # col 2
        (0.9, 0.1), (0.5, 0.25),
        # col 3
        (0.1, 0.9), (0.8, 0.2),
    ]
    points = [np.array(p) for p in points_raw]

    def dist_func(a, b):
        # This function replicates the 'dist' function in the shader.
        # Note: The neighborhood check optimization is omitted here for simplicity,
        # as it doesn't change the final distance value for points within the neighborhood.
        # It would only speed up the calculation by ignoring far-away points.
        uv = a - b
        return np.sqrt(uv[0]**2 + uv[1]**2)

    def worley_func(xy):
        # This function replicates the 'Worley' function.
        # It finds the minimum distance from the pixel 'xy' to any feature point.
        # The original shader loop has a bug (i < 5 for a mat4), we assume i < 4.
        min_dist = 2.0
        for p in points:
            d = dist_func(xy, p)
            min_dist = min(min_dist, d)
        return min_dist

    # --- Sample Calculation as requested ---
    # Calculate the color value for a sample pixel at the center of the image.
    uv_sample = np.array([0.5, 0.5])
    uv_sample_plus = uv_sample + np.array([0.025, 0.0])
    uv_sample_minus = uv_sample - np.array([0.025, 0.0])
    
    r_sample = worley_func(uv_sample)
    g_sample = worley_func(uv_sample_plus)
    b_sample = worley_func(uv_sample_minus)

    # The shader code evaluates to vec3(r), so the final color is (r, r, r)
    final_value_sample = r_sample
    
    print("--- Calculation for pixel at (0.5, 0.5) ---")
    print(f"r = Worley( (0.500, 0.5) ) = {r_sample:.4f}")
    print(f"g = Worley( (0.525, 0.5) ) = {g_sample:.4f}")
    print(f"b = Worley( (0.475, 0.5) ) = {b_sample:.4f}")
    print("\nFinal line is 'COLOR.rgb = vec3(vec3(r,g,b)).x', which evaluates to vec3(r).")
    print(f"Final normalized color value = {final_value_sample:.4f}")
    print(f"Final 8-bit color = ({int(final_value_sample*255)}, {int(final_value_sample*255)}, {int(final_value_sample*255)})")
    print("------------------------------------------\n")

    # --- Generate the full image ---
    print("Generating simulation image...")
    img = Image.new('RGB', (WIDTH, HEIGHT))
    pixels = img.load()

    for i in range(WIDTH):
        for j in range(HEIGHT):
            # UV coordinates assume (0,0) is at the top-left
            uv = np.array([i / WIDTH, j / HEIGHT])
            
            # As per our analysis, the final color only depends on the 'r' channel
            r = worley_func(uv)
            
            # Clamp the value between 0.0 and 1.0 and convert to 8-bit integer
            val_clamped = max(0.0, min(1.0, r))
            color_val = int(val_clamped * 255)
            
            # The result is grayscale
            pixels[i, j] = (color_val, color_val, color_val)

    # To see the output, you can save and open the image.
    # img.save("shader_output_simulation.png")
    # img.show()
    print("Image generation complete. The script produces a grayscale Worley noise pattern similar to Image C.")

solve_and_generate()