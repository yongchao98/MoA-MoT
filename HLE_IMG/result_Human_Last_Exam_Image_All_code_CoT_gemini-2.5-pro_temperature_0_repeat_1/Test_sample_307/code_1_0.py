import requests
from PIL import Image
import numpy as np
from io import BytesIO
import math

def solve_image_ratio():
    """
    Downloads, analyzes the image of Wigner 3-j symbols, and calculates
    the ratio of the maximum to the minimum infinity-norm.
    """
    try:
        # 1. Load the Image
        url = "https://i.imgur.com/9a5a2pA.png"
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for bad status codes
        img = Image.open(BytesIO(response.content)).convert('RGB')
        img_array = np.array(img)
    except requests.exceptions.RequestException as e:
        print(f"Error downloading the image: {e}")
        return

    # 2. Divide into Sub-images and Calculate Norms
    grid_size = 3
    img_width, img_height, _ = img_array.shape
    cell_width = img_width // grid_size
    cell_height = img_height // grid_size

    infinity_norms = []
    print("Calculating infinity-norm for each plot...")

    for i in range(grid_size):
        for j in range(grid_size):
            # Define the box for the current sub-image
            left = j * cell_width
            top = i * cell_height
            right = left + cell_width
            bottom = top + cell_height
            
            sub_img_array = img_array[top:bottom, left:right]
            
            max_magnitude = 0.0
            
            # 3. Iterate through pixels to find the max magnitude (infinity-norm)
            for pixel_row in sub_img_array:
                for pixel in pixel_row:
                    R, G, B = int(pixel[0]), int(pixel[1]), int(pixel[2])
                    
                    # Filter out black borders and white background/labels
                    # A simple sum check is effective for this image
                    if sum(pixel) < 60 or sum(pixel) > 720:
                        continue
                    
                    # 4. Define a Metric for Infinity-Norm
                    # Calculate Euclidean distance from white (255, 255, 255)
                    magnitude = math.sqrt((255 - R)**2 + (255 - G)**2 + (255 - B)**2)
                    
                    if magnitude > max_magnitude:
                        max_magnitude = magnitude
            
            plot_number = i * grid_size + j + 1
            print(f"Plot {plot_number}: Infinity-norm = {max_magnitude:.2f}")
            infinity_norms.append(max_magnitude)

    # 5. Find Max/Min and Ratio
    if not infinity_norms or min(infinity_norms) == 0:
        print("Could not determine norms or minimum norm is zero.")
        return
        
    max_norm = max(infinity_norms)
    min_norm = min(infinity_norms)
    ratio = max_norm / min_norm

    # 6. Calculate and Output the Ratio
    print("\n--- Calculation ---")
    print(f"Maximum norm found: {max_norm:.4f}")
    print(f"Minimum norm found: {min_norm:.4f}")
    print(f"Ratio of maximum to minimum norm is:")
    print(f"{max_norm:.4f} / {min_norm:.4f} = {ratio:.4f}")

solve_image_ratio()