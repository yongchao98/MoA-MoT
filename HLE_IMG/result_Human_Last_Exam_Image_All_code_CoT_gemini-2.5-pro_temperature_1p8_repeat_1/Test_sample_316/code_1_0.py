import numpy as np
import urllib.request
from PIL import Image
from scipy.signal import convolve2d
from itertools import permutations

def solve_ca_mapping():
    """
    Main function to solve the cellular automata mapping problem.
    """
    # --- 1. Load and Digitize All Patterns ---
    
    # URL of the image provided by the user.
    # In a real scenario, this would be a local file path 'image.png'.
    # For a self-contained script, we load it from the web.
    url = 'https://i.imgur.com/KPGYq6E.png'
    try:
        with urllib.request.urlopen(url) as img_url:
            img = Image.open(img_url).convert('L')
    except:
        print("Failed to download image. Please ensure you have an internet connection,")
        print("or place the image file as 'image.png' in the same directory.")
        return

    img_array = np.array(img)

    # Panel coordinates and sizes based on image layout (185x185, 45x45 panels, 1px border)
    panel_size = 45
    border_size = 1
    y_coords = [border_size + i * (panel_size + border_size) for i in range(4)]
    x_coords = [border_size + i * (panel_size + border_size) for i in range(4)]

    panel_labels = [
        ['A', 'B', 'C', 'D'], ['E', 'F', 'G', 'H'],
        ['1', '2', '3', '4'], ['5', '6', '7', '8']
    ]
    patterns = {}
    for r, row in enumerate(panel_labels):
        for c, label in enumerate(row):
            panel_img = img_array[y_coords[r]:y_coords[r]+panel_size, x_coords[c]:x_coords[c]+panel_size]
            panel_binary = (panel_img < 128).astype(int)
            if label.isalpha():
                patterns[label] = panel_binary[2:43, 2:43]  # Crop to 41x41
            else:
                patterns[label] = panel_binary[2:42, 2:42]  # Crop to 40x40

    # --- 2. Reverse-Engineer the Rules (A-H) ---
    def get_rule_map(rule_int):
        rule_bin = bin(rule_int)[2:].zfill(6)
        return np.array([int(c) for c in reversed(rule_bin)])

    def simulate_ca(rule_int):
        steps, size = 20, 81
        grid = np.zeros((size, size), dtype=int)
        center = size // 2
        grid[center, center] = 1
        kernel = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]])
        rule_map = get_rule_map(rule_int)
        for _ in range(steps):
            sums = convolve2d(grid, kernel, mode='same', boundary='fill', fillvalue=0)
            grid = rule_map[sums]
        return grid[center-steps:center+steps+1, center-steps:center+steps+1]

    ah_labels = "ABCDEFGH"
    found_rules = {}
    for rule_int in range(64):
        sim_pattern = simulate_ca(rule_int)
        for label in ah_labels:
            if label not in found_rules and np.array_equal(sim_pattern, patterns[label]):
                found_rules[label] = rule_int

    # --- 3. Establish the Mapping via Permutation Search ---
    labels_ah_sorted = sorted(found_rules.keys())
    rules_sorted = [found_rules[label] for label in labels_ah_sorted]

    # RAT[i, s] = result of rule i for sum s
    rule_application_table = np.array([get_rule_map(r) for r in rules_sorted])
    
    # The set of 6 possible 8-bit vectors resulting from applying all 8 rules
    valid_cols = {tuple(rule_application_table[:, s]) for s in range(6)}

    labels_18 = [str(i) for i in range(1, 9)]
    solution_perm = None

    for perm in permutations(labels_18):
        # Test hypothesis: Rule A -> perm[0], B -> perm[1], etc.
        is_valid_perm = True
        # Create a map from sorted alpha label index to permuted numeric label
        perm_map = {i: perm[i] for i in range(8)}
        
        for y in range(40):
            for x in range(40):
                # Observed 8-bit vector for pixel (y,x) based on current permutation
                obs_vec = tuple(patterns[perm_map[i]][y, x] for i in range(8))
                
                if obs_vec not in valid_cols:
                    is_valid_perm = False
                    break
            if not is_valid_perm:
                break
        
        if is_valid_perm:
            solution_perm = perm
            break
            
    # --- 4. Final Output ---
    final_mapping = {label: solution_perm[i] for i, label in enumerate(labels_ah_sorted)}
    
    result_list = [final_mapping[label] for label in "ABCDEFGH"]
    
    print("Rule mapping (A-H) to patterns (1-8):")
    for i, label in enumerate("ABCDEFGH"):
        print(f"Rule {label} corresponds to pattern {result_list[i]}")

    final_answer = "{" + ",".join(result_list) + "}"
    print("\nFinal answer in the required format:")
    print(final_answer)

solve_ca_mapping()