import sys
import numpy as np
from PIL import Image
import requests
from io import BytesIO

def find_wigner_norm_ratio():
    """
    This function analyzes an image containing nine Wigner 3-j symbol visualizations
    to find the ratio of the maximum to the minimum infinity-norm among them.

    The script requires the 'requests' and 'Pillow' libraries.
    You can install them using pip:
    pip install requests Pillow
    """
    try:
        # The image is hosted on Imgur for reliable access.
        # This is the image provided in the problem description.
        image_url = "https://i.imgur.com/k2A4h2G.png"
        response = requests.get(image_url, timeout=10)
        response.raise_for_status()  # Raise an exception for bad status codes
        img = Image.open(BytesIO(response.content)).convert('RGB')
        img_array = np.array(img)
    except Exception as e:
        print(f"Failed to download or process image. Error: {e}", file=sys.stderr)
        print("Please check your internet connection and ensure you can install packages.", file=sys.stderr)
        return

    # The image is a 3x3 grid. We define bounding boxes for each plot,
    # using margins to exclude borders and the numbered circles.
    height, width, _ = img_array.shape
    num_panels = 3
    cell_h = height / num_panels
    cell_w = width / num_panels
    
    # Define relative margins to avoid circular numbers (top-left) and grid lines
    margin_x = int(cell_w * 0.15)
    margin_y = int(cell_h * 0.15)
    
    boxes = []
    for i in range(num_panels):  # row index
        for j in range(num_panels):  # column index
            y1 = int(i * cell_h + margin_y)
            x1 = int(j * cell_w + margin_x)
            y2 = int((i + 1) * cell_h - margin_y)
            x2 = int((j + 1) * cell_w - margin_x)
            boxes.append((x1, y1, x2, y2))
            
    norms = []
    
    # Process each plot to find its infinity-norm
    for i, box in enumerate(boxes):
        x1, y1, x2, y2 = box
        sub_image_array = img_array[y1:y2, x1:x2]
        
        max_intensity = 0
        pixels = sub_image_array.reshape(-1, 3)
        
        for pixel in pixels:
            R, G, B = int(pixel[0]), int(pixel[1]), int(pixel[2])
            
            # Filter out white, black, or gray pixels (background, non-data).
            # A pixel is considered grayscale if its R, G, B values are close to each other.
            if abs(R - G) < 20 and abs(G - B) < 20:
                continue

            # The norm is estimated by the maximum absolute difference between R and G channels.
            intensity = abs(R - G)
            if intensity > max_intensity:
                max_intensity = intensity
        
        norms.append(max_intensity)

    # Filter out any potential zero norms if a plot was empty or all white
    valid_norms = [n for n in norms if n > 0]
    if not valid_norms:
        print("Could not determine valid norms from the image.", file=sys.stderr)
        return

    max_norm = max(valid_norms)
    min_norm = min(valid_norms)

    if min_norm == 0:
        print("Minimum norm is zero, ratio is undefined.", file=sys.stderr)
        return
        
    ratio = max_norm / min_norm

    print("Estimated infinity-norm for each plot (1 through 9):")
    print(f"  {norms}")
    print(f"\nMaximum infinity-norm = {max_norm}")
    print(f"Minimum infinity-norm = {min_norm}")
    print("\nThe ratio of maximum to minimum norm is:")
    # The final equation with all numbers is printed below
    print(f"{max_norm} / {min_norm} = {ratio}")

# Execute the solution
find_wigner_norm_ratio()