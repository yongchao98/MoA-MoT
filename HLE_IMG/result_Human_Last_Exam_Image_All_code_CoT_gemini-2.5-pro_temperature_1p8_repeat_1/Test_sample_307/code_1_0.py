import numpy as np
from PIL import Image
import os

def solve_wigner_ratio():
    """
    Analyzes the provided image of Wigner 3-j symbols to find the ratio
    of the maximum to the minimum infinity-norm.

    This script assumes the image is named 'image.png' and is in the same
    directory. It requires the Pillow and NumPy libraries. You can install them using:
    pip install Pillow numpy
    """
    image_filename = 'image.png'
    
    # Check if the image file exists
    if not os.path.exists(image_filename):
        print(f"Error: The image file '{image_filename}' was not found.")
        print("Please make sure the image is saved in the same directory as this script.")
        return

    # Load the image and convert it to a NumPy array for processing
    try:
        with Image.open(image_filename) as img:
            img_rgb = img.convert('RGB')
            data = np.array(img_rgb)
    except Exception as e:
        print(f"Error opening or processing the image: {e}")
        return

    # Define the pixel boundaries for each of the 9 subplots.
    # The image is 628x628, with 2px black borders.
    # The slices are defined as [start_row:end_row, start_col:end_col].
    y_slices = [slice(1, 209), slice(211, 419), slice(421, 628)]
    x_slices = [slice(1, 209), slice(211, 419), slice(421, 628)]

    infinity_norms = []
    
    # Iterate through each subplot region
    for y_slice in y_slices:
        for x_slice in x_slices:
            # Extract the subplot data
            subplot_data = data[y_slice, x_slice, :]
            
            # To prevent data type issues, convert to a floating-point type
            subplot_float = subplot_data.astype(np.float32)
            
            # Extract Red and Green channels
            r_channel = subplot_float[:, :, 0]
            g_channel = subplot_float[:, :, 1]
            
            # Calculate magnitude as abs(R - G) for each pixel
            magnitude_map = np.abs(r_channel - g_channel)
            
            # The infinity-norm is the maximum magnitude in the subplot
            norm = np.max(magnitude_map)
            infinity_norms.append(norm)

    # Find the maximum and minimum norms from the list
    max_norm = max(infinity_norms)
    min_norm = min(infinity_norms)

    # Calculate the ratio
    ratio = max_norm / min_norm

    # Print the detailed results
    print("Infinity-norm for each symbol visualization:")
    for i, norm in enumerate(infinity_norms):
        print(f"  Symbol {i + 1}: {norm:.2f}")

    print(f"\nMaximum norm found: {max_norm:.2f}")
    print(f"Minimum norm found: {min_norm:.2f}")

    # Print the final equation as requested
    print("\nRatio of maximum to minimum norm:")
    print(f"{max_norm:.2f} / {min_norm:.2f} = {ratio}")

if __name__ == '__main__':
    solve_wigner_ratio()
