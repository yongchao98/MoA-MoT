import numpy as np
import imageio.v2 as imageio

def solve_ca_mapping():
    """
    Solves the cellular automata mapping problem by correlating pattern densities.
    """
    try:
        # 1. Load the image from the URL
        url = "https://i.imgur.com/kS5v2oE.png"
        img_color = imageio.imread(url)
    except Exception as e:
        print(f"Error loading image: {e}")
        print("Please ensure you have an internet connection and the imageio and numpy libraries are installed (`pip install imageio numpy`).")
        return

    # 2. Process the image
    # Convert to grayscale if it's not already, and then to binary (black=1, white=0)
    if img_color.ndim == 3:
        # Use luminance formula for grayscale conversion before thresholding
        img_gray = np.dot(img_color[...,:3], [0.2989, 0.5870, 0.1140])
    else:
        img_gray = img_color
    
    img_binary = (img_gray < 128).astype(int)

    # 3. Segment and Calculate Densities
    H, W = img_binary.shape
    if H < 600 or W < 600:
        print("Image dimensions are smaller than expected. Cannot proceed.")
        return
        
    panel_H, panel_W = H // 4, W // 4
    
    panel_labels = [
        ['A', 'B', 'C', 'D'],
        ['E', 'F', 'G', 'H'],
        ['1', '2', '3', '4'],
        ['5', '6', '7', '8']
    ]

    densities = {}
    # Trim the border around each panel to focus on the content
    border_trim = 10 

    for i in range(4):
        for j in range(4):
            label = panel_labels[i][j]
            y_start, x_start = i * panel_H, j * panel_W
            y_end, x_end = y_start + panel_H, x_start + panel_W
            
            # Extract the panel and trim its borders
            panel = img_binary[y_start + border_trim : y_end - border_trim,
                               x_start + border_trim : x_end - border_trim]
            
            # The density is the sum of black pixels
            densities[label] = np.sum(panel)

    # Separate densities into the two groups
    labels_A_H = "ABCDEFGH"
    labels_1_8 = "12345678"
    
    densities_A_H = {k: densities[k] for k in labels_A_H}
    densities_1_8 = {k: densities[k] for k in labels_1_8}

    # 4. Correlate by sorting
    # Sort panel labels based on their calculated densities
    sorted_A_H = sorted(densities_A_H.keys(), key=lambda k: densities_A_H[k])
    sorted_1_8 = sorted(densities_1_8.keys(), key=lambda k: densities_1_8[k])

    # Create the mapping by pairing the sorted lists
    # e.g., lowest density A-H maps to lowest density 1-8
    mapping = {ah_label: num_label for ah_label, num_label in zip(sorted_A_H, sorted_1_8)}

    # 5. Format and Print the Output
    # Order the results according to the A-H sequence
    result_list = [mapping[label] for label in labels_A_H]
    
    # Construct the final answer string
    final_answer = "{" + ",".join(result_list) + "}"
    
    print(final_answer)

solve_ca_mapping()