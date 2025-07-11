import numpy as np
from PIL import Image
import requests
import io

def solve_ca_mapping():
    """
    Solves the cellular automata mapping problem by comparing pattern densities.
    """
    try:
        # Step 1: Fetch the image from the URL
        image_url = "https://i.imgur.com/kYRhR9L.png"
        response = requests.get(image_url)
        response.raise_for_status()  # Raise an exception for bad status codes
        img = Image.open(io.BytesIO(response.content))
    except requests.exceptions.RequestException as e:
        print(f"Error fetching image: {e}")
        return

    # Step 2: Define the bounding boxes for each sub-image (A-H, 1-8)
    # Format: (left, upper, right, lower)
    boxes = {
        'A': (3, 3, 164, 164), 'B': (170, 3, 331, 164), 'C': (337, 3, 498, 164), 'D': (504, 3, 665, 164),
        'E': (3, 170, 164, 331), 'F': (170, 170, 331, 331), 'G': (337, 170, 498, 331), 'H': (504, 170, 665, 331),
        '1': (3, 337, 164, 498), '2': (170, 337, 331, 498), '3': (337, 337, 498, 498), '4': (504, 337, 665, 498),
        '5': (3, 504, 164, 665), '6': (170, 504, 331, 665), '7': (337, 504, 498, 665), '8': (504, 504, 665, 665)
    }

    group1_labels = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    group2_labels = ['1', '2', '3', '4', '5', '6', '7', '8']

    densities_ah = []
    densities_18 = []

    def extract_grid(img_crop, grid_dim):
        """Extracts a CA grid by sampling pixels from an image crop."""
        img_crop_gray = img_crop.convert('L')
        sub_img_w, sub_img_h = img_crop_gray.size
        cell_w = sub_img_w / grid_dim
        cell_h = sub_img_h / grid_dim
        grid = np.zeros((grid_dim, grid_dim), dtype=int)
        for r in range(grid_dim):
            for c in range(grid_dim):
                px, py = int((c + 0.5) * cell_w), int((r + 0.5) * cell_h)
                pixel_val = img_crop_gray.getpixel((px, py))
                if pixel_val < 128:  # Black pixel
                    grid[r, c] = 1
        return grid

    # Step 3: Process Group 1 (A-H), calculate densities
    grid_dim_ah = 41
    for label in group1_labels:
        img_crop = img.crop(boxes[label])
        grid = extract_grid(img_crop, grid_dim_ah)
        density = np.mean(grid)
        densities_ah.append({'label': label, 'density': density})

    # Step 4: Process Group 2 (1-8), calculate densities
    grid_dim_18 = 40
    for label in group2_labels:
        img_crop = img.crop(boxes[label])
        grid = extract_grid(img_crop, grid_dim_18)
        density = np.mean(grid)
        densities_18.append({'label': label, 'density': density})

    # Step 5: Sort both groups by density to find the correspondence
    densities_ah.sort(key=lambda x: x['density'])
    densities_18.sort(key=lambda x: x['density'])

    # Step 6: Create the mapping and prepare the final result
    mapping = {}
    for i in range(len(densities_ah)):
        mapping[densities_ah[i]['label']] = densities_18[i]['label']

    print("Calculated densities and sorted order:\n")
    print("Group 1 (A-H) sorted by density:")
    for item in densities_ah:
        print(f"  Pattern {item['label']}: {item['density']:.4f}")

    print("\nGroup 2 (1-8) sorted by density:")
    for item in densities_18:
        print(f"  Pattern {item['label']}: {item['density']:.4f}")
    
    print("\nBased on the density ranking, the established mapping is:")
    result_list = []
    for label in group1_labels:
        mapped_label = mapping[label]
        result_list.append(mapped_label)
        print(f"  Rule {label} -> Pattern {mapped_label}")
    
    # Final Answer Formatting
    final_answer_str = ",".join(result_list)
    print("\nFinal Answer:")
    print(f"{{{final_answer_str}}}")

# Execute the solution
solve_ca_mapping()