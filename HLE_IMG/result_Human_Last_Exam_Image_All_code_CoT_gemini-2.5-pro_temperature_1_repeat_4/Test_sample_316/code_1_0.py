import numpy as np
from PIL import Image
from scipy.signal import convolve2d
from itertools import permutations
import requests
from io import BytesIO

def solve_ca_puzzle():
    """
    This script solves the cellular automata mapping puzzle by digitizing the images,
    reverse-engineering the rules, and solving a system of linear equations to
    find the correct mapping.
    """

    # --- Part 1: Image Processing ---
    def image_to_grids(img_url):
        """
        Loads an image from a URL and extracts the 16 40x40 CA grids.
        """
        try:
            response = requests.get(img_url)
            response.raise_for_status()
            img = Image.open(BytesIO(response.content)).convert('L')
        except Exception as e:
            print(f"Error: Could not load or process the image from URL. {e}")
            return None

        img_array = np.array(img)
        grids = {}
        labels = [c for c in 'ABCDEFGH'] + [str(i) for i in range(1, 9)]
        GRID_SIZE = 40
        SUBPLOT_SIZE = 250
        CROP_START, CROP_END = 24, 226
        content_size = CROP_END - CROP_START

        for i in range(4):
            for j in range(4):
                label = labels[i * 4 + j]
                top, left = i * SUBPLOT_SIZE, j * SUBPLOT_SIZE
                subplot = img_array[top:top + SUBPLOT_SIZE, left:left + SUBPLOT_SIZE]
                content = subplot[CROP_START:CROP_END, CROP_START:CROP_END]
                
                grid = np.zeros((GRID_SIZE, GRID_SIZE), dtype=int)
                cell_size = content_size / GRID_SIZE
                for r in range(GRID_SIZE):
                    for c in range(GRID_SIZE):
                        r_s, r_e = int(r * cell_size), int((r + 1) * cell_size)
                        c_s, c_e = int(c * cell_size), int((c + 1) * cell_size)
                        if np.mean(content[r_s:r_e, c_s:c_e]) < 128:
                            grid[r, c] = 1
                grids[label] = grid
        return grids

    # --- Part 2: Cellular Automata Simulation and Rule Identification ---
    def run_ca(rule_int, steps=20, grid_size=41):
        """
        Simulates a 5-neighbor totalistic CA from a single central cell.
        """
        rule_table = np.array([(rule_int >> i) & 1 for i in range(6)])
        grid = np.zeros((grid_size, grid_size), dtype=int)
        grid[grid_size // 2, grid_size // 2] = 1
        kernel = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], dtype=int)
        
        for _ in range(steps):
            sums = convolve2d(grid, kernel, mode='same', boundary='fill', fillvalue=0)
            grid = rule_table[sums]
        return grid

    def find_rules_for_A_H(target_grids):
        """
        Finds the rule for each pattern A-H by simulation and matching.
        """
        generated_patterns = {}
        sim_grid_size = 41
        center, crop_size = sim_grid_size // 2, 40
        
        for i in range(16): # 2^4 possibilities for b2,b3,b4,b5
            rule_int = (i << 2) + 2 # Rule format: b5 b4 b3 b2 1 0
            pattern = run_ca(rule_int, steps=20, grid_size=sim_grid_size)
            cropped_pattern = pattern[center-crop_size//2:center+crop_size//2, 
                                      center-crop_size//2:center+crop_size//2]
            generated_patterns[rule_int] = cropped_pattern

        rule_for_label = {}
        used_rules = set()
        for label in sorted(target_grids.keys()):
            target = target_grids[label]
            min_diff, best_rule = float('inf'), -1
            for rule_int, pattern in generated_patterns.items():
                if rule_int in used_rules: continue
                diff = np.sum(np.abs(target - pattern))
                if diff < min_diff:
                    min_diff, best_rule = diff, rule_int
            
            if best_rule != -1:
                rule_for_label[label] = best_rule
                used_rules.add(best_rule)
        return rule_for_label

    # --- Part 3: Solving the Mapping ---
    def solve_mapping(rule_for_label, grids_1_to_8):
        """
        Solves for the mapping by iterating through permutations and solving a linear system.
        """
        labels_A_H = sorted(rule_for_label.keys())
        rules_list = [rule_for_label[label] for label in labels_A_H]
        rule_matrix = np.array([[(r >> i) & 1 for i in range(6)] for r in rules_list])
        
        labels_1_8 = sorted(grids_1_to_8.keys(), key=int)
        observed_counts = np.array([np.sum(grids_1_to_8[label]) for label in labels_1_8])
        
        for p_indices in permutations(range(8)):
            permuted_rule_matrix = rule_matrix[list(p_indices), :]
            M_aug = np.vstack([permuted_rule_matrix, np.ones(6)])
            C_aug = np.append(observed_counts, 40*40)
            
            try:
                H, residuals, _, _ = np.linalg.lstsq(M_aug, C_aug, rcond=None)
            except np.linalg.LinAlgError:
                continue

            if residuals.size > 0 and residuals[0] < 1e-6 and np.all(H >= -1e-6) and np.allclose(H, np.round(H)):
                mapping = {}
                for i in range(8):
                    rule_label = labels_A_H[p_indices[i]]
                    count_label = str(i + 1)
                    mapping[rule_label] = count_label
                
                result_list = [mapping[label] for label in labels_A_H]
                return "{" + ",".join(result_list) + "}"
        return "Solution not found"

    # --- Main Execution ---
    image_url = "https://i.imgur.com/kLqgP2L.png"
    all_grids = image_to_grids(image_url)
    if all_grids is None: return

    grids_A_H = {k: v for k, v in all_grids.items() if k in 'ABCDEFGH'}
    rule_for_label = find_rules_for_A_H(grids_A_H)
    
    grids_1_to_8 = {k: v for k, v in all_grids.items() if k not in 'ABCDEFGH'}
    
    final_mapping = solve_mapping(rule_for_label, grids_1_to_8)
    print(final_mapping)

solve_ca_puzzle()