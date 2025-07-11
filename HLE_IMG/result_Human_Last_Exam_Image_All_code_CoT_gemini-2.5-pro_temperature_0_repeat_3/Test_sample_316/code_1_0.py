import numpy as np
from PIL import Image
from scipy.signal import convolve2d
from scipy.optimize import linear_sum_assignment
import itertools
import requests
import io

def solve_cellular_automata_mapping():
    """
    Solves the cellular automata mapping problem by identifying rules from
    long-term patterns and matching them to single-step patterns via density analysis.
    """

    # --- 1. Image Processing ---
    def load_and_crop_patterns():
        """
        Loads the image from a URL and crops the 16 patterns.
        """
        # The image is provided in the context, let's use its URL.
        # This makes the script runnable without manual downloads.
        image_url = "https://i.imgur.com/3a2tY2R.png"
        try:
            response = requests.get(image_url)
            response.raise_for_status()
            img = Image.open(io.BytesIO(response.content)).convert('L')
        except requests.exceptions.RequestException as e:
            print(f"Error: Could not download the image from {image_url}.")
            print(f"Please ensure you have an internet connection.")
            return None, None

        img_array = np.array(img)
        # Convert to binary: 1 for black, 0 for white
        binary_array = (img_array < 128).astype(int)

        patterns_A_H = {}
        patterns_1_8 = {}
        
        labels_A_H = [chr(ord('A') + i) for i in range(8)]
        labels_1_8 = [str(i + 1) for i in range(8)]

        # Crop patterns A-H (41x41) from the top two rows
        for i in range(8):
            row, col = i // 4, i % 4
            label = labels_A_H[i]
            cell_size = 100
            cell = binary_array[row*cell_size:(row+1)*cell_size, col*cell_size:(col+1)*cell_size]
            center_r, center_c = cell_size // 2, cell_size // 2
            radius = 41 // 2
            r_start, r_end = center_r - radius, center_r + radius + 1
            c_start, c_end = center_c - radius, center_c + radius + 1
            patterns_A_H[label] = cell[r_start:r_end, c_start:c_end]

        # Crop patterns 1-8 (40x40) from the bottom two rows
        for i in range(8):
            row, col = (i // 4) + 2, i % 4
            label = labels_1_8[i]
            cell_size = 100
            cell = binary_array[row*cell_size:(row+1)*cell_size, col*cell_size:(col+1)*cell_size]
            center_r, center_c = cell_size // 2, cell_size // 2
            radius = 40 // 2
            r_start, r_end = center_r - radius, center_r + radius
            c_start, c_end = center_c - radius, center_c + radius
            patterns_1_8[label] = cell[r_start:r_end, c_start:c_end]
            
        return patterns_A_H, patterns_1_8

    # --- 2. Rule Identification ---
    def simulate_ca(rule_tuple, steps, size):
        """Simulates a 5-neighbor totalistic CA from a single central cell."""
        grid = np.zeros((size, size), dtype=int)
        center = size // 2
        grid[center, center] = 1
        
        kernel = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], dtype=int)
        rule_array = np.array(rule_tuple, dtype=int)

        for _ in range(steps):
            sums = convolve2d(grid, kernel, mode='same', boundary='fill', fillvalue=0)
            grid = rule_array[sums]
            
        return grid

    def find_rules_for_A_H(patterns_A_H):
        """Finds the rule for each pattern A-H by simulating candidate rules."""
        letter_to_rule = {}
        # From problem analysis, all rules have b0=0 and b1=1.
        # We search the 16 rules where b2, b3, b4, b5 vary.
        for b_tuple in itertools.product([0, 1], repeat=4):
            b2, b3, b4, b5 = b_tuple
            rule = (0, 1, b2, b3, b4, b5)
            
            sim_result = simulate_ca(rule, steps=20, size=41)
            
            for label, pattern in patterns_A_H.items():
                if np.array_equal(sim_result, pattern):
                    letter_to_rule[label] = rule
                    break
    
        return letter_to_rule

    # --- 3. Density Calculation ---
    def calculate_theoretical_density(rule_tuple, p_initial=0.5):
        """Calculates the expected output density for a given rule and initial density."""
        b0, b1, b2, b3, b4, b5 = rule_tuple
        p = p_initial
        q = 1 - p
        
        # p' = sum_{k=0 to 5} [ C(5, k) * p^k * q^(5-k) * b_k ]
        p_prime = (
            1 * (p**0)*(q**5) * b0 + 5 * (p**1)*(q**4) * b1 +
            10 * (p**2)*(q**3) * b2 + 10 * (p**3)*(q**2) * b3 +
            5 * (p**4)*(q**1) * b4 + 1 * (p**5)*(q**0) * b5
        )
        return p_prime

    # --- Main Execution Logic ---
    patterns_A_H, patterns_1_8 = load_and_crop_patterns()
    if not patterns_A_H:
        return

    letter_to_rule = find_rules_for_A_H(patterns_A_H)
    
    if len(letter_to_rule) != 8:
        print("Error: Could not identify all 8 rules from patterns A-H.")
        print("This might be due to an issue with the simulation or image processing.")
        return

    measured_densities = {label: np.mean(p) for label, p in patterns_1_8.items()}
    
    # --- 4. Optimal Assignment ---
    letters = sorted(letter_to_rule.keys())
    numbers = sorted(measured_densities.keys(), key=int)
    
    theo_densities = {l: calculate_theoretical_density(letter_to_rule[l]) for l in letters}

    cost_matrix = np.zeros((8, 8))
    for i, letter in enumerate(letters):
        for j, number in enumerate(numbers):
            cost_matrix[i, j] = abs(theo_densities[letter] - measured_densities[number])
            
    row_ind, col_ind = linear_sum_assignment(cost_matrix)
    
    final_mapping = {}
    for i in range(8):
        letter = letters[row_ind[i]]
        number = numbers[col_ind[i]]
        final_mapping[letter] = number
        
    # --- 5. Report the Result ---
    result_list = [final_mapping[chr(ord('A') + i)] for i in range(8)]
    result_str = "{" + ",".join(result_list) + "}"
    print("Found the following mapping for (A, B, C, D, E, F, G, H):")
    print(result_str)

solve_cellular_automata_mapping()