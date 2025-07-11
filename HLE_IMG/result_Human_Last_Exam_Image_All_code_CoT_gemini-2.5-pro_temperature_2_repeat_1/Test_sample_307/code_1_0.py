import numpy as np
from PIL import Image
import requests
import io

def find_norm_ratio():
    """
    Analyzes an image of nine Wigner 3-j symbols to find the ratio
    of the maximum to the minimum infinity-norm.
    """
    try:
        # The image is hosted on Imgur for stable access.
        image_url = "https://i.imgur.com/T0bHnGe.png"
        response = requests.get(image_url)
        response.raise_for_status()  # Ensure the request was successful
        img = Image.open(io.BytesIO(response.content)).convert("RGB")
        img_array = np.array(img)
    except requests.exceptions.RequestException as e:
        print(f"Error fetching image: {e}")
        return

    # Define the bounding boxes for each of the 9 subplots.
    # The format is (x_start, y_start, x_end, y_end).
    # These coordinates were determined by inspecting the grid lines in the image.
    boxes = [
        # Row 1
        (5, 30, 183, 210), (185, 30, 365, 210), (367, 30, 545, 210),
        # Row 2
        (5, 212, 183, 389), (185, 212, 365, 389), (367, 212, 545, 389),
        # Row 3
        (5, 391, 183, 545), (185, 391, 365, 545), (367, 391, 545, 545),
    ]

    infinity_norms = []
    white_pixel = np.array([255, 255, 255])

    for i, (x1, y1, x2, y2) in enumerate(boxes):
        # Extract the subplot from the main image array
        subplot_array = img_array[y1:y2, x1:x2]
        
        # Calculate the Euclidean distance of each pixel from pure white.
        # This distance serves as a proxy for the absolute value of the Wigner symbol.
        distances = np.linalg.norm(subplot_array - white_pixel, axis=2)
        
        # The infinity-norm is the maximum distance found in the subplot.
        max_dist = np.max(distances)
        infinity_norms.append(max_dist)

    # Find the overall maximum and minimum norms from the list
    max_norm = max(infinity_norms)
    min_norm = min(infinity_norms)

    # Calculate the final ratio
    ratio = max_norm / min_norm

    print("This script calculates the ratio of the maximum to minimum infinity-norm.")
    print("The infinity-norm for each heatmap is estimated by finding the maximum color distance from white.\n")
    print(f"The maximum norm is found in subplot #{infinity_norms.index(max_norm) + 1}, with a value of {max_norm:.2f}")
    print(f"The minimum norm is found in subplot #{infinity_norms.index(min_norm) + 1}, with a value of {min_norm:.2f}\n")
    print("The final calculation is:")
    print(f"{max_norm:.2f} / {min_norm:.2f} = {ratio}")

if __name__ == "__main__":
    find_norm_ratio()