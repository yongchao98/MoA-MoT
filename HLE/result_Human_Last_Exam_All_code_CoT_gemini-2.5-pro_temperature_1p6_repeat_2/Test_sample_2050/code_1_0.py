import numpy as np

def calculate_ground_state_energy():
    """
    Calculates and explains the ground state energy for a system of 4 photons
    in 7 coupled nonlinear cavities in the hard-core limit (U -> infinity).
    """
    # System parameters
    N = 7  # Number of cavities
    M = 4  # Number of photons

    # The problem is equivalent to finding the ground state of M non-interacting
    # fermions on a ring of N sites. The single-particle energy levels for the
    # hopping part of the Hamiltonian are given by E_q = -2J * cos(2*pi*q/N),
    # where q is an integer from 0 to N-1.

    # To find the M-particle ground state energy, we fill the M lowest-energy
    # single-particle states.

    # Calculate energies for all q states.
    # Energies are proportional to -cos(2*pi*q/N), so we find the largest cos values.
    q_vals = np.arange(N)
    cos_k_vals = np.cos(2 * np.pi * q_vals / N)

    # Get the indices that would sort cos(k) in descending order (i.e., ascending energy)
    sorted_indices = np.argsort(cos_k_vals)[::-1]

    # The q values of the M lowest-energy states are the first M in the sorted list.
    occupied_q = q_vals[sorted_indices[:M]]

    # The total energy is E_gs = E_onsite + E_interaction + E_hopping.
    # E_onsite = sum(omega * n_i) = omega * M = 4*omega.
    # E_interaction = 0 in the U -> infinity limit for the ground state.
    # E_hopping = sum over occupied states of E_q = sum_{q_occ} -2J * cos(2*pi*q_occ/N).

    print("The system is equivalent to 4 non-interacting fermions on a 7-site ring.")
    print(f"The ground state is found by filling the {M} lowest single-particle energy states.")
    print(f"The q-values for the {M} lowest energy states are: {list(occupied_q)} (where q=6 corresponds to q=-1, etc.).")
    
    # We construct the expression for the ground state energy E_gs.
    print("\nThe full expression for the ground state energy E_gs is:")
    
    # Group terms for a cleaner expression, e.g., 2*cos(...)
    term_counts = {}
    for q in occupied_q:
        # Use the q representation with the smallest absolute value, e.g., N-1 becomes -1.
        key = q if q <= N/2 else q - N
        # Group degenerate terms by the absolute value of q
        abs_key = abs(key)
        term_counts[abs_key] = term_counts.get(abs_key, 0) + 1

    # Build the string for the cosine sum part of the equation
    cos_terms_str = []
    # Sort keys to print in a canonical order (0, 1, 2, ...)
    sorted_keys = sorted(term_counts.keys())
    for key in sorted_keys:
        count = term_counts[key]
        
        # Format the term: "cos(...)", "1", "2*cos(...)"
        if key == 0:
            term_str = "1" # cos(0)=1
        else:
            term_str = f"cos(2*pi*{key}/{N})"
        
        if count > 1:
            term_str = f"{count}*{term_str}"
        
        cos_terms_str.append(term_str)
    
    # Assemble the final equation string
    # E_gs = M*omega - 2*J * ( sum of cos terms )
    hopping_energy_expr = " + ".join(cos_terms_str)
    final_expression = f"E_gs = {M}*omega - 2*J * ({hopping_energy_expr})"
    
    # Print each part of the equation as requested
    print(final_expression)
    
    # Calculate the numerical value of the hopping energy for verification.
    cos_sum = np.sum(cos_k_vals[occupied_q])
    hopping_energy_factor = -2 * cos_sum
    
    print("\nNumerically, the ground state energy is approximately:")
    print(f"E_gs = {M}*omega + ({hopping_energy_factor:.6f})*J")


# Run the calculation and print the results
calculate_ground_state_energy()
