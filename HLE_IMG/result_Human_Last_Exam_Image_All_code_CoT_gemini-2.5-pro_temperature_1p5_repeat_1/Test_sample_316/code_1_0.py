import numpy as np
from PIL import Image
import os

def solve_ca_mapping():
    """
    Solves the cellular automata mapping problem by calculating and matching pattern densities.
    """
    image_filename = 'image.png'

    # Check if the image file exists
    if not os.path.exists(image_filename):
        print(f"Error: The image file '{image_filename}' was not found.")
        print("Please download the image from the problem description and save it as 'image.png' in the same directory as this script.")
        return

    # Load the image and convert to grayscale
    img = Image.open(image_filename).convert('L')
    img_array = np.array(img)

    # Define labels and panel grid coordinates (column, row)
    labels_g1 = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    labels_g2 = [str(i) for i in range(1, 9)]
    
    panel_coords = {
        'A': (0, 0), 'B': (1, 0), 'C': (2, 0), 'D': (3, 0),
        'E': (0, 1), 'F': (1, 1), 'G': (2, 1), 'H': (3, 1),
        '1': (0, 2), '2': (1, 2), '3': (2, 2), '4': (3, 2),
        '5': (0, 3), '6': (1, 3), '7': (2, 3), '8': (3, 3),
    }

    # Image and panel dimensions
    img_w, img_h = img.size
    panel_w = img_w / 4
    panel_h = img_h / 4

    # --- Step 1: Calculate densities for Group 1 (A-H) ---
    densities_g1 = []
    print("Analyzing Group 1 (A-H)...")
    for label in labels_g1:
        c, r = panel_coords[label]
        # Center of the panel
        cx = int((c + 0.5) * panel_w)
        cy = int((r + 0.5) * panel_h)
        
        # Crop a 41x41 box for patterns A-H
        size = 41
        half = size // 2
        box = (cx - half, cy - half, cx + half + 1, cy + half + 1)
        
        # Extract the pattern and convert to binary (1 for black, 0 for white)
        cropped_array = img_array[box[1]:box[3], box[0]:box[2]]
        binary_array = (cropped_array < 128).astype(int)
        density = np.mean(binary_array)
        densities_g1.append({'label': label, 'density': density})
        print(f"  Pattern {label}: Density = {density:.4f}")

    # --- Step 2: Calculate densities for Group 2 (1-8) ---
    densities_g2 = []
    print("\nAnalyzing Group 2 (1-8)...")
    for label in labels_g2:
        c, r = panel_coords[label]
        cx = int((c + 0.5) * panel_w)
        cy = int((r + 0.5) * panel_h)
        
        # Crop a 40x40 box for patterns 1-8
        size = 40
        half = size // 2
        box = (cx - half, cy - half, cx + half, cy + half)
        
        cropped_array = img_array[box[1]:box[3], box[0]:box[2]]
        binary_array = (cropped_array < 128).astype(int)
        density = np.mean(binary_array)
        densities_g2.append({'label': label, 'density': density})
        print(f"  Pattern {label}: Density = {density:.4f}")

    # --- Step 3: Sort both groups by density to find the mapping ---
    sorted_g1 = sorted(densities_g1, key=lambda x: x['density'])
    sorted_g2 = sorted(densities_g2, key=lambda x: x['density'])

    print("\n--- Density Rankings ---")
    print("Group 1 (low to high):", [item['label'] for item in sorted_g1])
    print("Group 2 (low to high):", [item['label'] for item in sorted_g2])

    # Create the mapping based on sorted lists
    mapping = {}
    for i in range(len(sorted_g1)):
        label_g1 = sorted_g1[i]['label']
        label_g2 = sorted_g2[i]['label']
        mapping[label_g1] = label_g2
        
    # --- Step 4: Report the final result ---
    result_order = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    final_map_numbers = [mapping[label] for label in result_order]

    print("\n--- Final Mapping ---")
    for label in result_order:
        print(f"Rule {label} -> Pattern {mapping[label]}")

    print("\nFinal Answer Format:")
    final_answer_string = "{" + ",".join(final_map_numbers) + "}"
    print(final_answer_string)


if __name__ == '__main__':
    solve_ca_mapping()