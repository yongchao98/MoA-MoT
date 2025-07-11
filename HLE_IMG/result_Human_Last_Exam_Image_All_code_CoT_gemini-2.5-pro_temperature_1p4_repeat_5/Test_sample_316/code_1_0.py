import numpy as np
from PIL import Image

def solve_ca_mapping():
    """
    Solves the cellular automata mapping problem by correlating pattern densities.
    """
    try:
        # Load the image and convert it to a binary numpy array (black=1, white=0)
        img = Image.open('image.png').convert('L')
        img_array = np.array(img)
        # A threshold of 128 is used to binarize the image, handling any anti-aliasing.
        bin_array = (img_array < 128).astype(int)
    except FileNotFoundError:
        print("Error: 'image.png' not found.")
        print("Please ensure the image file is in the same directory as the script.")
        return

    # Image dimensions and panel size calculation
    img_h, img_w = bin_array.shape
    x_step = img_w / 4.0
    y_step = img_h / 4.0

    def get_pixel_count(row, col):
        """
        Crops a panel from the grid and counts the black pixels, avoiding borders and labels.
        """
        x_start = int(col * x_step)
        y_start = int(row * y_step)
        x_end = int((col + 1) * x_step)
        y_end = int((row + 1) * y_step)

        # Define margins to crop out labels and grid lines.
        # These values are tuned based on visual inspection of the image layout.
        y_margin_top = 38
        y_margin_bottom = 10
        
        top = y_start + y_margin_top
        bottom = y_end - y_margin_bottom
        
        # Panels A-H have labels on the top-right; Panels 1-8 on the top-left.
        # Adjust horizontal margins accordingly.
        if row < 2:  # Panels A-H
            x_margin_left = 10
            x_margin_right = 38
        else:  # Panels 1-8
            x_margin_left = 38
            x_margin_right = 10

        left = x_start + x_margin_left
        right = x_end - x_margin_right
        
        panel_data = bin_array[top:bottom, left:right]
        return np.sum(panel_data)

    panel_labels_top = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    panel_labels_bottom = ['1', '2', '3', '4', '5', '6', '7', '8']

    # Calculate black pixel counts for each panel
    counts_top = {label: get_pixel_count(i // 4, i % 4) for i, label in enumerate(panel_labels_top)}
    counts_bottom = {label: get_pixel_count(2 + i // 4, i % 4) for i, label in enumerate(panel_labels_bottom)}

    # Sort both groups by pixel count to find the correlation
    sorted_top = sorted(counts_top.items(), key=lambda item: item[1])
    sorted_bottom = sorted(counts_bottom.items(), key=lambda item: item[1])

    print("--- Analysis Results ---")
    print("\nPatterns A-H sorted by black pixel count (density):")
    for label, count in sorted_top:
        print(f"  Pattern {label}: {count} pixels")
    
    print("\nPatterns 1-8 sorted by black pixel count (density):")
    for label, count in sorted_bottom:
        print(f"  Pattern {label}: {count} pixels")

    # Create the mapping based on the sorted order
    mapping = {sorted_top[i][0]: sorted_bottom[i][0] for i in range(8)}
    
    print("\nEstablished Mapping (A-H -> 1-8):")
    for label_A_H, label_1_8 in mapping.items():
        print(f"  Rule {label_A_H} -> Pattern {label_1_8}")

    # Prepare the final answer in the specified format
    result_list = [mapping[label] for label in panel_labels_top]
    result_str = ",".join(result_list)
    
    print("\n--- Final Answer ---")
    print(f"The final mapping for A,B,C,D,E,F,G,H is:")
    print(f"{{{result_str}}}")

# Execute the solution
solve_ca_mapping()