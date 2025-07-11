import numpy as np
from PIL import Image

def solve_wigner_ratio():
    """
    Calculates the ratio of the maximum to the minimum infinity-norm
    from an image of nine Wigner 3-j symbol heatmaps.
    """
    try:
        # Load the image. Ensure it's in the same directory as the script.
        # The image should be saved as 'wigner_symbols.png'
        img = Image.open('wigner_symbols.png').convert('RGB')
    except FileNotFoundError:
        print("Error: 'wigner_symbols.png' not found.")
        print("Please save the image from the problem description as 'wigner_symbols.png' in the correct directory.")
        return

    pixels = np.array(img)
    
    # Image is a 3x3 grid. Calculate the size of each subplot.
    box_height, box_width = pixels.shape[0] // 3, pixels.shape[1] // 3

    norms = []
    # Loop through the 3x3 grid of subplots
    for i in range(3):
        for j in range(3):
            # Define the bounding box for the current subplot
            y_start, y_end = i * box_height, (i + 1) * box_height
            x_start, x_end = j * box_width, (j + 1) * box_width

            subplot_pixels = pixels[y_start:y_end, x_start:x_end]

            max_intensity = 0
            
            # Iterate through pixels in the subplot to find the max intensity
            for y in range(subplot_pixels.shape[0]):
                for x in range(subplot_pixels.shape[1]):
                    r, g, b = subplot_pixels[y, x]
                    
                    # Filter out grayscale pixels (background, borders, numbers)
                    # A pixel is considered grayscale if R, G, and B values are close to each other.
                    if abs(int(r) - int(g)) < 15 and abs(int(g) - int(b)) < 15:
                        continue
                        
                    # The intensity (a proxy for the norm) is the difference between red and green channels.
                    intensity = abs(int(r) - int(g))
                    
                    if intensity > max_intensity:
                        max_intensity = intensity
            
            norms.append(max_intensity)

    # Find the maximum and minimum norms from the list
    max_norm = max(norms)
    min_norm = min(norms)

    if min_norm == 0:
        print("Could not calculate ratio: the minimum norm is zero.")
    else:
        # Calculate the final ratio
        ratio = max_norm / min_norm
        
        print(f"Maximum infinity-norm found: {max_norm}")
        print(f"Minimum infinity-norm found: {min_norm}")
        print("\nThe ratio is calculated as follows:")
        # The final equation is printed as requested
        print(f"{max_norm} / {min_norm} = {ratio}")

if __name__ == '__main__':
    solve_wigner_ratio()