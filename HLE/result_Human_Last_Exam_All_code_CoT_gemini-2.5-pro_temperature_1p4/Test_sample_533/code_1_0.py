import itertools
from collections import defaultdict

def solve_h2_fock_space():
    """
    Constructs and decomposes the Fock space for H2 in a minimal basis.
    """
    # 1. Define the minimal basis of spin-orbitals
    # (spatial_symmetry, spin)
    spin_orbitals = [
        ('g', 'up'),    # sigma_g(alpha)
        ('g', 'down'),  # sigma_g(beta)
        ('u', 'up'),    # sigma_u(alpha)
        ('u', 'down')   # sigma_u(beta)
    ]
    num_spin_orbitals = len(spin_orbitals)

    print("Step 1: Define the basis")
    print(f"Minimal basis consists of {num_spin_orbitals} spin-orbitals:")
    for so in spin_orbitals:
        print(f"  - |σ_{so[0]}, {so[1]}⟩")
    print("-" * 30)

    # 2. Construct the full Fock space (all 2^4 = 16 states)
    # Each state is an occupancy vector, e.g., (1, 1, 0, 0) means
    # the first two spin-orbitals are occupied.
    fock_space_states = list(itertools.product([0, 1], repeat=num_spin_orbitals))

    print("Step 2: Construct the Fock space")
    print(f"The Fock space contains 2^{num_spin_orbitals} = {len(fock_space_states)} total electronic states.")
    print("-" * 30)

    # 3. Group states by conserved quantum numbers (N, Parity)
    # The dictionary will store M_S values for each (N, Parity) group.
    grouped_states = defaultdict(list)

    for state in fock_space_states:
        # Calculate N (number of electrons)
        n_electrons = sum(state)

        # Calculate M_S (total spin projection)
        m_s = 0.0
        for i, occupied in enumerate(state):
            if occupied:
                m_s += 0.5 if spin_orbitals[i][1] == 'up' else -0.5

        # Calculate g/u parity (product of parities of occupied orbitals)
        # 'g' parity is +1, 'u' parity is -1
        parity_val = 1
        for i, occupied in enumerate(state):
            if occupied:
                if spin_orbitals[i][0] == 'u':
                    parity_val *= -1
        parity = 'g' if parity_val == 1 else 'u'
        
        # In a minimal sigma-only basis, all states are Sigma_plus.
        # We only need to track N, S, and g/u parity.

        grouped_states[(n_electrons, parity)].append(m_s)

    print("Step 3: Group states by N and Parity")
    print("The 16 states are grouped by (Number of Electrons, Parity).")
    print("The lists show the M_S values found in each group:")
    for (n, p), ms_list in sorted(grouped_states.items()):
        print(f"  - N={n}, Parity='{p}': M_s values = {sorted(ms_list)}")
    print("-" * 30)

    # 4. Identify unique Hilbert spaces by deducing S
    hilbert_spaces = set()
    space_counts_by_N = defaultdict(int)

    for (n, parity), ms_values in sorted(grouped_states.items()):
        # Create a sorted list of M_S values to peel off multiplets
        ms_list = sorted(list(ms_values), reverse=True)
        while ms_list:
            # The highest M_S value in the list corresponds to S
            s = ms_list[0]
            
            # This (N, S, Parity) tuple uniquely identifies a Hilbert space
            space_term = f"^{int(2*s)+1}Σ_{parity}"
            # Adding a '+' superscript as all are of this type.
            if int(2*s)+1 > 0: # Avoid superscript for vacuum
                space_term += "⁺"
            
            space_label = (n, s, parity, space_term)
            hilbert_spaces.add(space_label)
            space_counts_by_N[n] += 1
            
            # Remove the M_S components of this newly found multiplet
            for k in range(int(2 * s) + 1):
                m = s - k
                try:
                    ms_list.remove(m)
                except ValueError:
                    # This should not happen in a correct physical model
                    print(f"Error: Incomplete multiplet found for S={s} in group N={n}, P='{parity}'")

    print("Step 4: Identify unique symmetry-adapted Hilbert spaces")
    print("By analyzing the M_S values, we identify the following spaces (N, S, Parity):")
    sorted_spaces = sorted(list(hilbert_spaces))
    for n, s, parity, term in sorted_spaces:
        # A small fix for the vacuum state term symbol
        if n == 0: term = "¹Σg⁺ (Vacuum)"
        print(f"  - N={n}, S={s}, Parity='{parity}'  => corresponds to space of type {term}")
    print("-" * 30)

    # 5. Count the total number of spaces
    total_spaces = len(hilbert_spaces)
    print("Step 5: Count the maximum number of symmetry-adapted Hilbert spaces")
    print("The Fock space decomposes into the direct sum of these spaces.")
    
    count_strings = [str(c) for n, c in sorted(space_counts_by_N.items())]
    print(f"Number of spaces for N=0, 1, 2, 3, 4 are: {', '.join(count_strings)}")
    print(f"The final equation for the total number of spaces is:")
    print(f"Total Spaces = {' + '.join(count_strings)} = {total_spaces}")
    print("-" * 30)
    print(f"The maximum number of symmetry-adapted Hilbert spaces is {total_spaces}.")
    
    return total_spaces

# Run the analysis
final_answer = solve_h2_fock_space()
print(f"<<<{final_answer}>>>")