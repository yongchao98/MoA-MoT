import math
from PIL import Image
import requests
from io import BytesIO

def solve():
    """
    Finds the ratio of the maximum to the minimum infinity-norm among the nine
    Wigner 3-j symbols visualized in the image.
    """
    # Step 1: Load the image.
    # We will try to download the image from a known URL to ensure we have the exact file.
    # If it fails, it will look for a local file named 'image.png'.
    image_path = 'image.png'
    try:
        # This is the original image from the prompt context.
        url = "https://i.stack.imgur.com/G5lYw.png"
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for bad status codes
        img = Image.open(BytesIO(response.content))
    except Exception:
        try:
            # Fallback to a local file if the download fails
            img = Image.open(image_path)
        except FileNotFoundError:
            print(f"Error: Could not find or download '{image_path}'.")
            print("Please make sure the image file is in the same directory as the script.")
            return

    # --- Helper Functions ---

    def get_color_magnitude(rgb):
        """Calculates the Euclidean distance of an RGB color from white."""
        r, g, b = rgb[:3]
        return math.sqrt((255 - r)**2 + (255 - g)**2 + (255 - b)**2)

    def is_data_pixel(rgb):
        """
        Determines if a pixel is part of the data. This filters out the
        black frame, white background, and grayish numbered circles.
        A pixel is considered data if it's sufficiently colored (not grayscale).
        """
        r, g, b = rgb[:3]
        # Ignore pure white pixels (background)
        if r > 250 and g > 250 and b > 250:
            return False
        # Ignore dark pixels (black frame and numbers)
        if r < 30 and g < 30 and b < 30:
            return False
        # Ignore grayscale pixels (e.g., the numbered circles) by checking saturation
        if max(r, g, b) - min(r, g, b) < 15:
            return False
        return True

    # --- Main Analysis ---

    width, height = img.size
    cell_w = width // 3
    cell_h = height // 3
    plot_norms = []

    # Step 2: Iterate through each of the 9 plots.
    for i in range(3):  # Row index
        for j in range(3):  # Column index
            # Define the bounding box for the current plot cell
            left, upper = j * cell_w, i * cell_h
            right, lower = (j + 1) * cell_w, (i + 1) * cell_h
            cell_img = img.crop((left, upper, right, lower))

            max_magnitude_in_plot = 0
            # Step 3: Find the max magnitude (infinity-norm) in the current plot.
            for pixel in cell_img.getdata():
                if is_data_pixel(pixel):
                    magnitude = get_color_magnitude(pixel)
                    if magnitude > max_magnitude_in_plot:
                        max_magnitude_in_plot = magnitude
            
            plot_norms.append(max_magnitude_in_plot)

    # Step 4: Find the overall maximum and minimum norms from the list.
    if not plot_norms or min(p for p in plot_norms if p > 0) == 0:
        print("Could not find valid data in the image.")
        return
        
    max_norm = max(plot_norms)
    # Ensure min_norm is not zero
    min_norm = min(p for p in plot_norms if p > 0)

    # For context, identify which plots had the min and max norms
    max_plot_index = plot_norms.index(max_norm) + 1
    min_plot_index = plot_norms.index(min_norm) + 1

    print("--- Analysis Results ---")
    print(f"Maximum norm (most intense color) found in plot {max_plot_index}: {max_norm:.4f}")
    print(f"Minimum norm (least intense color) found in plot {min_plot_index}: {min_norm:.4f}")
    
    # Step 5: Calculate and print the ratio.
    ratio = max_norm / min_norm
    print("\n--- Final Calculation ---")
    print(f"Ratio = Maximum Norm / Minimum Norm")
    print(f"Ratio = {max_norm:.4f} / {min_norm:.4f} = {ratio:.4f}")

if __name__ == '__main__':
    solve()