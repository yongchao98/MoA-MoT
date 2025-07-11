import numpy as np
from PIL import Image
import requests
from io import BytesIO

def solve():
    """
    This script solves the cellular automata mapping problem by following these steps:
    1.  Downloads and processes the problem's image file.
    2.  Extracts the eight patterns from Group 2 (labeled 1-8) and calculates their density (number of black pixels).
    3.  Sorts these eight images based on their densities in descending order.
    4.  Uses a pre-determined productivity ranking of the rules A-H. This ranking was established by identifying the specific rule for each pattern A-H, analyzing their theoretical productivity, and using the visual density of the A-H patterns themselves to resolve any ties.
    5.  Maps the sorted list of rules to the sorted list of images to establish the final correspondence.
    6.  Formats and prints the final answer.
    """
    
    # URL of the image file for the problem
    image_url = "https://i.imgur.com/b9t3kUQ.png"
    try:
        response = requests.get(image_url)
        response.raise_for_status()  # Ensure the download was successful
        # Open the image and convert to grayscale
        img = Image.open(BytesIO(response.content)).convert('L')
    except Exception as e:
        print(f"Fatal Error: Could not download or process the image from {image_url}.")
        print(f"Reason: {e}")
        return

    # Define geometric properties of the image grid based on visual inspection
    cell_size = 107      # Size of each of the 16 grid cells
    pattern_offset = 14  # Inner margin to the top-left of the 80x80 pattern
    pattern_size = 80    # Pixel size of the pattern area in the image
    ca_size = 40         # Dimension of the cellular automaton grid (40x40)

    def extract_pattern(row, col):
        """
        Extracts a single CA pattern from the main image grid.
        The function crops the 80x80 pixel area and downsamples it to a 40x40 numpy array
        of 0s (white) and 1s (black).
        """
        # Calculate the top-left corner of the 80x80 crop box
        left = col * cell_size + pattern_offset
        top = row * cell_size + pattern_offset
        
        # Crop the image
        crop_box = (left, top, left + pattern_size, top + pattern_size)
        pattern_image = img.crop(crop_box)
        
        # Binarize and downsample the cropped image into a numpy array
        pattern_array = np.zeros((ca_size, ca_size), dtype=np.int32)
        for r in range(ca_size):
            for c in range(ca_size):
                # Sample the top-left pixel of each 2x2 block
                pixel_value = pattern_image.getpixel((c * 2, r * 2))
                # Threshold to determine black (1) or white (0)
                if pixel_value < 128:
                    pattern_array[r, c] = 1
        return pattern_array

    # Calculate the density for each image in Group 2 (1-8)
    image_densities = {}
    for i in range(8):
        # Images 1-4 are in row 2, 5-8 are in row 3 (0-indexed)
        row = 2 + (i // 4)
        col = i % 4
        image_label = str(i + 1)
        pattern = extract_pattern(row, col)
        image_densities[image_label] = np.sum(pattern)

    # Sort the image labels '1' through '8' by their density in descending order
    sorted_image_labels = sorted(image_densities, key=image_densities.get, reverse=True)

    # This is the crucial ranking of rules from most productive to least productive.
    # It was determined through simulation, analysis, and a tie-breaking heuristic.
    rule_productivity_order = ['H', 'F', 'E', 'B', 'G', 'D', 'C', 'A']

    # Map the most productive rule to the densest image, and so on.
    final_mapping = {}
    for i in range(len(rule_productivity_order)):
        rule_label = rule_productivity_order[i]
        image_label = sorted_image_labels[i]
        final_mapping[rule_label] = image_label
    
    # Construct the final result string in the specified A,B,C... order
    result_list = [final_mapping[letter] for letter in "ABCDEFGH"]
    
    # Print the answer in the required format
    print(f"{{{','.join(result_list)}}}")

# Run the solver function
solve()