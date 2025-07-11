import numpy as np
from PIL import Image

def get_panels(image_path='image.png'):
    """
    Crops the main image into 16 panels (A-H, 1-8).
    Each panel is resized to a 41x41 grid and converted to a binary matrix
    where Black=1 and White=0.
    """
    try:
        img = Image.open(image_path).convert('L')
    except FileNotFoundError:
        print(f"Error: The image file was not found at '{image_path}'.")
        print("Please download the image from the prompt and save it as 'image.png' in the same folder as this script.")
        return None

    # Automatically find the content area by trimming the white border
    img_array = np.array(img)
    non_white_pixels = np.where(img_array < 250)
    if non_white_pixels[0].size > 0:
        y_min, y_max = np.min(non_white_pixels[0]), np.max(non_white_pixels[0])
        x_min, x_max = np.min(non_white_pixels[1]), np.max(non_white_pixels[1])
        img = img.crop((x_min, y_min, x_max, y_max))

    width, height = img.size
    panel_w, panel_h = width // 4, height // 4

    labels = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H',
              '1', '2', '3', '4', '5', '6', '7', '8']
    panels = {}
    
    for i, label in enumerate(labels):
        row, col = i // 4, i % 4
        left, top = col * panel_w, row * panel_h
        
        # Crop each panel, adding an inner margin to remove grid lines
        margin = int(panel_w * 0.05)
        panel_img = img.crop((left + margin, top + margin, left + panel_w - margin, top + panel_h - margin))
        
        # Resize to 41x41, appropriate for t=20 evolution
        resized_panel = panel_img.resize((41, 41), Image.Resampling.LANCZOS)
        panel_array = np.array(resized_panel)
        
        # Convert to binary matrix (thresholding grayscale values)
        panels[label] = (panel_array < 128).astype(int)
        
    return panels

def run_ca(rule_int, steps, size):
    """
    Simulates a 5-neighbor totalistic CA for a given rule integer.
    The simulation starts from a single central cell.
    """
    grid = np.zeros((size, size), dtype=int)
    grid[size // 2, size // 2] = 1
    
    # A rule is defined by the outcome for each possible sum (0 to 5)
    rule_map = {(i): (rule_int >> i) & 1 for i in range(6)}

    for _ in range(steps):
        # Pad the grid to handle boundaries easily
        padded_grid = np.pad(grid, 1, mode='constant')
        new_grid = np.zeros_like(grid)
        for r in range(size):
            for c in range(size):
                # Calculate the sum over the 5-cell Von Neumann neighborhood
                s = (padded_grid[r + 1, c + 1] +  # Cell itself
                     padded_grid[r, c + 1]     +  # Top
                     padded_grid[r + 2, c + 1] +  # Bottom
                     padded_grid[r + 1, c]     +  # Left
                     padded_grid[r + 1, c + 2])   # Right
                
                new_grid[r, c] = rule_map.get(s, 0)
        grid = new_grid
    return grid

def solve():
    """
    Main function to execute the plan:
    1. Digitize patterns from the image.
    2. Find rules for A-H by simulating and matching.
    3. Score rules and count pixels in patterns 1-8.
    4. Match the sorted lists to find the final mapping.
    """
    print("Step 1: Digitizing patterns from image...")
    all_patterns = get_panels('image.png')
    if all_patterns is None:
        return

    # Cache all 32 possible simulation results (for rules where b0=0)
    simulation_cache = {}
    for i in range(32):
        rule_int = i * 2  # b0 is 0, so rule is even. b1-b5 vary.
        simulation_cache[rule_int] = run_ca(rule_int, steps=20, size=41)

    print("Step 2: Identifying rules for patterns A-H by finding the best match...")
    group1_labels = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    found_rules = {}
    for label in group1_labels:
        target_grid = all_patterns[label]
        min_diff = float('inf')
        best_rule = -1
        for rule_int, sim_grid in simulation_cache.items():
            diff = np.sum(np.abs(sim_grid - target_grid)) # L1 distance
            if diff < min_diff:
                min_diff = diff
                best_rule = rule_int
        found_rules[label] = best_rule

    # Use binomial coefficients as weights for the density score
    weights = {1: 5, 2: 10, 3: 10, 4: 5, 5: 1} # C(5,k) for k=1..5
    rule_scores = {}
    for label, rule_int in found_rules.items():
        score = sum(weights[k] * ((rule_int >> k) & 1) for k in range(1, 6))
        rule_scores[label] = score

    print("Step 3: Calculating densities...")
    group2_labels = ['1', '2', '3', '4', '5', '6', '7', '8']
    black_counts = {label: np.sum(all_patterns[label]) for label in group2_labels}

    # Match by sorting scores and counts
    sorted_group1 = sorted(group1_labels, key=lambda label: rule_scores[label])
    sorted_group2 = sorted(group2_labels, key=lambda label: black_counts[label])

    mapping = {g1: g2 for g1, g2 in zip(sorted_group1, sorted_group2)}

    # Construct the final answer in the required format
    final_answer_list = [mapping[label] for label in group1_labels]
    result_str = ",".join(final_answer_list)
    final_output_format = f"{{{result_str}}}"

    print("\n--- Summary ---")
    print("Rules (A-H) sorted by predicted density score:")
    print("  ", ", ".join([f"{l}({rule_scores[l]})" for l in sorted_group1]))
    print("Patterns (1-8) sorted by actual black pixel count:")
    print("  ", ", ".join([f"{l}({black_counts[l]})" for l in sorted_group2]))

    print("\n--- Final Mapping ---")
    for i, label in enumerate(group1_labels):
        print(f"Rule {label} corresponds to Pattern {final_answer_list[i]}")

    print("\nFinal Answer in specified format:")
    print(final_output_format)
    print(f"<<<{final_output_format}>>>")

if __name__ == '__main__':
    solve()