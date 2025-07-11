import numpy as np
import matplotlib.pyplot as plt

# Define the feature points from the shader
POINTS_MAT = np.array([
    [0.1, 0.1, 0.5, 0.5],
    [0.8, 0.8, 0.25, 0.5],
    [0.9, 0.1, 0.5, 0.25],
    [0.1, 0.9, 0.8, 0.2]
])

# Unpack points into a list of [x, y] coordinates
points_list = []
for i in range(4):
    points_list.append(np.array([POINTS_MAT[i, 0], POINTS_MAT[i, 1]]))
    points_list.append(np.array([POINTS_MAT[i, 2], POINTS_MAT[i, 3]]))

def chebyshev_dist(p1, p2):
    """
    Calculates Chebyshev distance (L∞ norm) to produce squares/rectangles.
    This replaces the original Euclidean distance function.
    """
    return np.max(np.abs(p1 - p2))

def worley(xy, points):
    """
    Finds the minimum Chebyshev distance from a point xy to any of the feature points.
    """
    min_dist = 10.0
    for p in points:
        min_dist = min(min_dist, chebyshev_dist(xy, p))
    return min_dist

def fragment_for_image_A(uv, points):
    """
    Simulates the logic required to generate Image A.
    """
    # 1. Calculate Worley distance for R, G, B channels at offset coordinates
    #    to create chromatic aberration.
    dist_r = worley(uv, points)
    dist_g = worley(uv + np.array([0.01, 0]), points) # Offset for G channel
    dist_b = worley(uv - np.array([0.01, 0]), points) # Offset for B channel

    # 2. Apply a threshold to create binary shapes (black on white).
    #    A pixel is black (0) if the distance is less than the threshold,
    #    and white (1) otherwise.
    threshold = 0.05
    r = 1.0 if dist_r > threshold else 0.0
    g = 1.0 if dist_g > threshold else 0.0
    b = 1.0 if dist_b > threshold else 0.0
    
    # The colors in the image (Cyan, Magenta, Yellow) are subtractive colors,
    # which appear when you overlay colored filters. On a screen (additive color),
    # a black shape for the red channel on a white background corresponds to
    # showing only Cyan light.
    # R channel shape -> (0, 1, 1) -> Cyan
    # G channel shape -> (1, 0, 1) -> Magenta
    # B channel shape -> (1, 1, 0) -> Yellow
    # All overlap -> (0, 0, 0) -> Black
    return np.array([r, g, b])

# Generate the image
width, height = 512, 512
image = np.zeros((height, width, 3))

for y in range(height):
    for x in range(width):
        # Convert pixel coordinates to UV coordinates (0 to 1)
        uv = np.array([x / width, y / height])
        image[y, x] = fragment_for_image_A(uv, points_list)

# The problem asks to print the equation. We will print the final color
# for a sample point, e.g., the center of the image (0.5, 0.5), which is
# near one of the feature points.
uv_center = np.array([0.5, 0.5])
color_at_center = fragment_for_image_A(uv_center, points_list)

print("To generate image A, the shader logic must be modified.")
print("The key changes are using Chebyshev distance, enabling chromatic aberration, and thresholding the result.")
print(f"For example, at the central point UV=(0.5, 0.5):")
dist_r = worley(uv_center, points_list)
dist_g = worley(uv_center + np.array([0.01, 0]), points_list)
dist_b = worley(uv_center - np.array([0.01, 0]), points_list)
print(f"  R distance = {dist_r:.4f}, G distance = {dist_g:.4f}, B distance = {dist_b:.4f}")
print(f"  With a threshold of 0.05, the final RGB color value is ({color_at_center[0]}, {color_at_center[1]}, {color_at_center[2]}), which is black.")

# To display the generated image (optional, for user verification)
# plt.imshow(image)
# plt.title("Python simulation for Image A")
# plt.show()
