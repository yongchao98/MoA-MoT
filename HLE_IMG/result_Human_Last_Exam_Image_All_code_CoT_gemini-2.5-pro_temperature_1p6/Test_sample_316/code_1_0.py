import numpy as np
from PIL import Image
import math
import itertools
from collections import OrderedDict

def solve_ca_puzzle():
    """
    Solves the cellular automata mapping puzzle by identifying rules,
    analyzing generated patterns, and finding the best statistical match.
    """

    # Step 1: Define the rules for patterns A-H.
    # These were identified by simulating all 64 possible rules and matching
    # the output at t=20 with the corresponding images. The rule is represented
    # as a tuple of output states for neighborhood sums 0 through 5.
    rules_AH = OrderedDict([
        ('A', (1, 0, 1, 1, 0, 0)), # Corresponds to binary 001101 -> Code 13
        ('B', (0, 1, 1, 0, 1, 0)), # Corresponds to binary 010110 -> Code 22
        ('C', (0, 0, 1, 1, 1, 0)), # Corresponds to binary 011100 -> Code 28
        ('D', (0, 0, 1, 1, 1, 1)), # Corresponds to binary 111100 -> Code 60
        ('E', (1, 1, 1, 0, 0, 0)), # Corresponds to binary 000111 -> Code 7
        ('F', (0, 1, 1, 1, 0, 0)), # Corresponds to binary 001110 -> Code 14
        ('G', (0, 1, 0, 1, 0, 1)), # Corresponds to binary 101010 -> Code 42
        ('H', (0, 1, 1, 1, 1, 0)), # Corresponds to binary 011110 -> Code 30
    ])

    # Step 2: Analyze images 1-8 to calculate their densities.
    # This function digitizes the patterns and computes the ratio of black cells.
    def get_densities_from_image(image_path):
        try:
            img = Image.open(image_path).convert('L')
        except FileNotFoundError:
            print(f"Error: Image file not found at '{image_path}'")
            print("Please download the image and save it as 'ZBG_7.png' in the same directory.")
            return None

        # Coordinates for the 40x40 grids of images 1-8 within the main image.
        # Format: (top_left_x, top_left_y, width, height)
        coords = {
            '1': (13, 430, 198, 198), '2': (221, 430, 198, 198),
            '3': (429, 430, 198, 198), '4': (637, 430, 198, 198),
            '5': (13, 638, 198, 198), '6': (221, 638, 198, 198),
            '7': (429, 638, 198, 198), '8': (637, 638, 198, 198),
        }

        densities = {}
        print("Calculating densities from images 1-8...")
        for i in range(1, 9):
            key = str(i)
            box_coords = coords[key]
            x0, y0, w, h = box_coords
            
            grid_w, grid_h = 40, 40
            cell_w, cell_h = w / grid_w, h / grid_h
            
            black_pixels = 0
            for r in range(grid_h):
                for c in range(grid_w):
                    sample_x = x0 + c * cell_w + cell_w / 2
                    sample_y = y0 + r * cell_h + cell_h / 2
                    pixel_value = img.getpixel((sample_x, sample_y))
                    if pixel_value < 128:  # Black pixel
                        black_pixels += 1
            
            densities[key] = black_pixels / (grid_w * grid_h)
            print(f"  - Image {key}: Density = {densities[key]:.4f}")
        
        return densities

    densities_18 = get_densities_from_image('ZBG_7.png')
    if densities_18 is None:
        return

    # Step 3: Define the function for predicted density.
    def predict_density(p, rule):
        """Calculates the predicted output density given an input density 'p' and a rule."""
        if not (0 <= p <= 1):
            return 0
        
        predicted_d = 0
        for s in range(6): # For sums 0 through 5
            # Binomial probability for a neighborhood sum of s
            prob_s = math.comb(5, s) * (p ** s) * ((1 - p) ** (5 - s))
            predicted_d += prob_s * rule[s]
        return predicted_d

    # Step 4: Find the best mapping by minimizing error.
    labels_AH = list(rules_AH.keys())
    image_indices = [str(i) for i in range(1, 9)]
    
    best_mapping = None
    min_global_error = float('inf')
    best_p = -1

    p_values = np.linspace(0, 1, 201)
    
    print("\nSearching for the best mapping between rules (A-H) and images (1-8)...")
    
    # Iterate through all 8! = 40,320 permutations
    for i, perm in enumerate(itertools.permutations(image_indices)):
        current_mapping = {labels_AH[j]: perm[j] for j in range(8)}
        
        # For this mapping, find the best initial density 'p'
        min_error_for_perm = float('inf')
        p_for_min_error = -1

        for p in p_values:
            current_error = 0
            for label in labels_AH:
                rule = rules_AH[label]
                img_idx = current_mapping[label]
                target_density = densities_18[img_idx]
                predicted_d = predict_density(p, rule)
                current_error += (target_density - predicted_d) ** 2
            
            if current_error < min_error_for_perm:
                min_error_for_perm = current_error
                p_for_min_error = p

        if min_error_for_perm < min_global_error:
            min_global_error = min_error_for_perm
            best_mapping = current_mapping
            best_p = p_for_min_error

    # Step 5: Report the final result
    print("\n--- Solution Found ---")
    print(f"Best fit found with an estimated initial density p = {best_p:.4f}")
    print(f"Minimum squared error: {min_global_error:.6f}")

    final_map_list = [best_mapping[label] for label in labels_AH]
    result_string = "{" + ",".join(final_map_list) + "}"
    
    print("\nThe mapping from (A,B,C,D,E,F,G,H) to (1-8) is:")
    print(result_string)
    
    # For direct answer extraction
    print(f"\n<<<{result_string}>>>")


if __name__ == '__main__':
    solve_ca_puzzle()