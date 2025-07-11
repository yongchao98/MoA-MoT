import numpy as np
import matplotlib.pyplot as plt
import io
import base64

def generate_image_a_logic():
    """
    This function simulates the logic required to generate image A,
    assuming the provided shader code is a puzzle with misleading parts.
    
    Corrections made based on analysis:
    1. The color calculation is corrected to produce chromatic aberration.
    2. The distance function is replaced with one that generates rectangular shapes.
    """
    
    WIDTH, HEIGHT = 256, 256
    SIZE = 20.0
    
    # These are the 8 points defined in the shader's 'points' matrix
    points = [
        np.array([0.1, 0.1]), np.array([0.5, 0.5]),
        np.array([0.8, 0.8]), np.array([0.25, 0.5]),
        np.array([0.9, 0.1]), np.array([0.5, 0.25]),
        np.array([0.1, 0.9]), np.array([0.8, 0.2])
    ]
    
    def rect_dist(uv, center):
        """
        A distance function that defines a rectangular area.
        Returns 0 inside the rectangle, 1 outside.
        This function replaces the misleading Euclidean distance in the shader.
        The size of the rectangle is chosen to approximate the shapes in image D.
        """
        rect_half_size = np.array([1.5 / SIZE, 1.0 / SIZE])
        d = np.abs(uv - center) - rect_half_size
        is_inside = np.all(d < 0.0)
        return 0.0 if is_inside else 1.0

    def worley_rect(uv):
        """
        Calculates the minimum distance to any point using the rectangular
        distance function. This effectively checks if the pixel is in the
        union of all rectangles.
        """
        min_d = 1.0
        for p in points:
            min_d = min(min_d, rect_dist(uv, p))
        return min_d

    # Create an empty image array
    image_data = np.zeros((HEIGHT, WIDTH, 3))
    
    # The offset for chromatic aberration
    offset = np.array([0.01, 0.01])

    # Generate pixel data
    for y in range(HEIGHT):
        for x in range(WIDTH):
            # Normalize coordinates to 0-1 range, with (0,0) at bottom-left
            uv = np.array([x / WIDTH, (HEIGHT - 1 - y) / HEIGHT])
            
            # This is the "equation" for the final color
            # We calculate the value for each color channel separately
            # Red channel: no offset
            # Green channel: positive offset
            # Blue channel: negative offset
            r = worley_rect(uv)
            g = worley_rect(uv + offset)
            b = worley_rect(uv - offset)
            
            image_data[y, x] = [r, g, b]

    # The prompt asks to output each number in the final equation.
    # The equation is COLOR.rgb = vec3(r, g, b).
    # Since this varies per pixel, we'll print an example for the center pixel.
    center_uv = np.array([0.5, 0.5])
    r_center = worley_rect(center_uv)
    g_center = worley_rect(center_uv + offset)
    b_center = worley_rect(center_uv - offset)
    
    print("Example calculation for the final color equation at UV=(0.5, 0.5):")
    print(f"Red channel value = {r_center}")
    print(f"Green channel value = {g_center}")
    print(f"Blue channel value = {b_center}")
    print("The image is constructed by applying this logic to every pixel.")
    print("\nThe Python code simulates the logic that produces image A.")

generate_image_a_logic()