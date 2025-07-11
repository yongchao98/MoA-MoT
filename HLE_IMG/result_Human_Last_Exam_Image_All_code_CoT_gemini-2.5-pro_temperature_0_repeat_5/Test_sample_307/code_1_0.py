import numpy as np
from PIL import Image
import urllib.request
from io import BytesIO

def solve_wigner_ratio():
    """
    This function calculates the ratio of the maximum to the minimum infinity-norm
    among the nine Wigner 3-j symbols visualized in the provided image.
    """
    # --- Step 1: Load the image ---
    # The script will attempt to download the image from a URL.
    # If you have it saved locally, you can replace the loading part with:
    # image_path = 'wigner_symbols.png'
    # img = Image.open(image_path).convert('RGB')
    try:
        url = "https://i.imgur.com/9a5a2pA.png"
        with urllib.request.urlopen(url) as response:
            img_data = response.read()
        img = Image.open(BytesIO(img_data)).convert('RGB')
        img_array = np.array(img, dtype=np.int32)
        print("Image loaded successfully.\n")
    except Exception as e:
        print(f"Error: Could not load the image. Please check the URL or save the image locally.")
        print(f"Details: {e}")
        return

    # --- Step 2: Define bounding boxes for the 9 sub-images ---
    # The image is 648x648, so each cell in the 3x3 grid is 216x216.
    # We'll use a margin to exclude the black grid lines.
    cell_size = 216
    margin = 10
    boxes = []
    for i in range(3):  # row index
        for j in range(3):  # column index
            y_start = i * cell_size + margin
            y_end = (i + 1) * cell_size - margin
            x_start = j * cell_size + margin
            x_end = (j + 1) * cell_size - margin
            boxes.append((y_start, y_end, x_start, x_end))

    # --- Step 3: Function to calculate the infinity-norm for a sub-image ---
    def get_infinity_norm(sub_image_array):
        # Filter out pure white background pixels.
        # A pixel is considered background if its R, G, and B values are all very high.
        # We create a mask to select only the colored pixels.
        is_not_background = np.sum(sub_image_array, axis=2) < 750
        colored_pixels = sub_image_array[is_not_background]

        if colored_pixels.shape[0] == 0:
            return 0  # No data pixels found in this sub-image

        # The value is modeled as (Green channel - Red channel).
        # The infinity-norm is the maximum of the absolute values.
        r_channel = colored_pixels[:, 0]
        g_channel = colored_pixels[:, 1]
        values = g_channel - r_channel
        norm = np.max(np.abs(values))
        return norm

    # --- Step 4: Calculate the norm for each of the 9 sub-images ---
    norms = []
    for i, box in enumerate(boxes):
        y_start, y_end, x_start, x_end = box
        sub_image = img_array[y_start:y_end, x_start:x_end]
        norm = get_infinity_norm(sub_image)
        norms.append(norm)

    print(f"Calculated infinity-norms for images 1 through 9:")
    print(f"{[int(n) for n in norms]}\n")

    # --- Step 5: Find the maximum and minimum norms ---
    max_norm = max(norms)
    min_norm = min(norms)

    print(f"Maximum norm found: {int(max_norm)}")
    print(f"Minimum norm found: {int(min_norm)}\n")

    # --- Step 6: Calculate and print the final ratio ---
    if min_norm == 0:
        print("Error: Minimum norm is zero, so the ratio is undefined.")
    else:
        ratio = max_norm / min_norm
        print("Final Calculation:")
        print(f"Ratio = Maximum Norm / Minimum Norm")
        print(f"Ratio = {int(max_norm)} / {int(min_norm)} = {ratio:.4f}")

if __name__ == '__main__':
    solve_wigner_ratio()