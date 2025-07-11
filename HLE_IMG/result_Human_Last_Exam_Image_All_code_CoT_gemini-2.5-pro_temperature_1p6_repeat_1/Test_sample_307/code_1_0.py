import numpy as np
from PIL import Image
import requests
import io

def solve_wigner_ratio():
    """
    This function performs the entire task:
    1. Downloads the image of Wigner 3-j symbols.
    2. Extracts the 9 matrices from the image.
    3. Calculates the infinity-norm for each matrix.
    4. Finds the ratio of the maximum norm to the minimum norm.
    5. Prints the results.
    """
    
    # URL of the image containing the Wigner 3-j symbols
    image_url = "https://i.imgur.com/B94j3oB.png"
    
    # --- Step 1: Download and open the image ---
    try:
        response = requests.get(image_url)
        response.raise_for_status() # Raise an exception for bad status codes
        img_data = response.content
        main_image = Image.open(io.BytesIO(img_data)).convert('RGB')
    except requests.exceptions.RequestException as e:
        print(f"Error: Could not download or open the image from {image_url}. {e}")
        return

    # --- Step 2: Extract matrices from the image grid ---
    matrices = []
    
    # Image grid parameters derived from visual inspection
    panel_size = 198
    gap = 7
    grid_offset = 6
    grid_step = panel_size + gap
    matrix_dim = 11
    block_size = panel_size / matrix_dim

    for panel_row in range(3):
        for panel_col in range(3):
            # Calculate the bounding box for the current panel
            left = grid_offset + panel_col * grid_step
            top = grid_offset + panel_row * grid_step
            right = left + panel_size
            bottom = top + panel_size
            
            panel_image = main_image.crop((left, top, right, bottom))
            
            # Initialize the matrix for the current panel
            matrix = np.zeros((matrix_dim, matrix_dim))
            
            for r in range(matrix_dim): # Matrix row
                for c in range(matrix_dim): # Matrix column
                    # Sample the color from the center of the block
                    px = c * block_size + block_size / 2
                    py = r * block_size + block_size / 2
                    
                    R, G, B = panel_image.getpixel((px, py))
                    
                    # Convert RGB color to a numerical value assuming a red-white-green colormap.
                    # The scaling is relative, which is sufficient for a ratio calculation.
                    value = (G - R) / 255.0
                    matrix[r, c] = value
            
            matrices.append(matrix)

    # --- Step 3: Calculate the infinity-norm for each matrix ---
    infinity_norms = []
    for i, matrix in enumerate(matrices):
        # The infinity-norm is the maximum absolute row sum.
        norm = np.max(np.sum(np.abs(matrix), axis=1))
        infinity_norms.append(norm)

    # --- Step 4: Find the max/min norms and calculate the ratio ---
    # Filter out any potential zero norms, although none are expected here.
    non_zero_norms = [norm for norm in infinity_norms if norm > 1e-9]
    
    if not non_zero_norms:
        print("Could not find any non-zero norms.")
        return

    max_norm = max(non_zero_norms)
    min_norm = min(non_zero_norms)
    ratio = max_norm / min_norm
    
    # --- Step 5: Print the final result as requested ---
    print("This script finds the ratio of the maximum to the minimum infinity-norm among the nine visualized Wigner 3-j symbols.")
    print("\nThe final calculation is:")
    print(f"{max_norm} / {min_norm} = {ratio}")

if __name__ == '__main__':
    solve_wigner_ratio()
