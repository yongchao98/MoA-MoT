import numpy as np
from PIL import Image
import requests
from io import BytesIO

def solve_ca_mapping():
    """
    Solves the cellular automata mapping problem by matching the density of black cells.
    """
    
    # 1. Download and load the image from the provided URL.
    image_url = "https://i.imgur.com/L7rP1bQ.png"
    try:
        response = requests.get(image_url)
        response.raise_for_status()  # Raise an exception for bad status codes
        # Open image and convert to grayscale for simplicity
        img = Image.open(BytesIO(response.content)).convert('L')
    except Exception as e:
        print(f"Error: Failed to download or process the image. {e}")
        return

    # 2. Define parameters for grid extraction.
    # The full image is 800x800, containing a 4x4 array of sub-images.
    SUB_IMG_SIZE = 200  # Each sub-image is 200x200 pixels
    GRID_DIM = 40       # The CA grid within each sub-image is 40x40
    CELL_SIZE = SUB_IMG_SIZE // GRID_DIM  # This means each CA cell is 5x5 pixels

    # 3. Define labels and coordinates for all 16 sub-images.
    labels_group1 = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    labels_group2 = [str(i) for i in range(1, 9)]
    all_labels = labels_group1 + labels_group2
    
    coords = {}
    for i in range(4):  # row index of sub-image
        for j in range(4):  # column index of sub-image
            label = all_labels[i * 4 + j]
            # Top-left corner of the sub-image
            coords[label] = (j * SUB_IMG_SIZE, i * SUB_IMG_SIZE)

    def get_black_cell_count(label):
        """
        Extracts the 40x40 grid for a given label and returns the count of black cells.
        """
        grid = np.zeros((GRID_DIM, GRID_DIM), dtype=int)
        x_offset, y_offset = coords[label]
        for r in range(GRID_DIM):
            for c in range(GRID_DIM):
                # Sample the pixel at the center of the 5x5 pixel cell
                px = x_offset + c * CELL_SIZE + (CELL_SIZE // 2)
                py = y_offset + r * CELL_SIZE + (CELL_SIZE // 2)
                
                pixel_value = img.getpixel((px, py))
                
                # In grayscale, black pixels have low values. We use 128 as a threshold.
                if pixel_value < 128:
                    grid[r, c] = 1
        return np.sum(grid)

    # 4. Calculate black cell counts for both groups.
    counts1 = {label: get_black_cell_count(label) for label in labels_group1}
    counts2 = {label: get_black_cell_count(label) for label in labels_group2}

    # 5. Sort each group's items by their counts to find the correspondence.
    sorted_group1 = sorted(counts1.items(), key=lambda item: item[1])
    sorted_group2 = sorted(counts2.items(), key=lambda item: item[1])

    # 6. Establish the mapping and print the underlying logic as requested.
    print("Mapping based on sorting by black cell counts:")
    mapping = {}
    for i in range(len(sorted_group1)):
        label1, count1 = sorted_group1[i]
        label2, count2 = sorted_group2[i]
        mapping[label1] = label2
        print(f"Rule {label1} (count: {count1:3d}) maps to Result {label2} (count: {count2:3d})")

    # 7. Format the final answer in the specified order (A, B, C, ..., H).
    result_list = [mapping[label] for label in labels_group1]
    result_str = "{" + ",".join(result_list) + "}"
    
    print("\nFinal Answer:")
    print(result_str)

solve_ca_mapping()