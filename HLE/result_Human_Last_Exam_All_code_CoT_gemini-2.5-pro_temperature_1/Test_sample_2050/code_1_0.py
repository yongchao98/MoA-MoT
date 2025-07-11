import numpy as np

def solve_hamiltonian():
    """
    Calculates the ground state energy for a system of 7 nonlinear optical cavities.

    The method involves:
    1. Simplifying the Hamiltonian in the U -> infinity limit, which corresponds to
       hard-core bosons (no two photons on the same site).
    2. Mapping the hard-core boson problem to a system of non-interacting spinless fermions.
    3. Finding the single-particle energy levels of the tight-binding Hamiltonian on a 7-site ring.
    4. Filling the lowest 4 energy levels to find the ground state energy of the hopping term.
    5. Combining this with the constant energy terms to get the total ground state energy.
    """
    # Number of cavities
    N = 7
    # Number of photons
    N_ph = 4

    # The Hamiltonian is H = H_const + H_int + H_hop
    # H_const = sum_i omega * a_i^dag a_i = omega * N_ph = 4*omega
    # H_int = sum_i U/2 * n_i(n_i-1). For U -> infinity, the ground state must have n_i in {0, 1}.
    # This makes H_int = 0. This is the hard-core boson condition.
    # H_hop = -J * sum_i (a_i^dag a_{i+1} + a_i^dag a_{i-1})
    # This system is equivalent to 4 spinless fermions on a 7-site ring.

    # The single-particle energies for the hopping Hamiltonian are E_m = -2*J*cos(2*pi*m/N)
    # where m = 0, +/-1, +/-2, ...
    # For N=7, the quantum numbers for the 7 distinct energy levels can be taken as m = 0, +/-1, +/-2, +/-3.

    # We need to fill the lowest 4 energy levels. The ground state is formed by occupying
    # the states corresponding to m = 0, +1, -1, and +2 (or -2, due to degeneracy).
    # E_hop = E_0 + E_{+1} + E_{-1} + E_{+2}
    # E_hop = -2*J * (cos(0) + 2*cos(2*pi/7) + cos(4*pi/7))

    # Let's calculate the numerical values.
    cos0 = np.cos(0)
    cos1 = np.cos(2 * np.pi / N)
    cos2 = np.cos(4 * np.pi / N)

    # The total coefficient for the -2*J term
    sum_of_cosines = cos0 + 2 * cos1 + cos2
    
    # The coefficient for the J term in the energy
    j_coefficient = -2 * sum_of_cosines

    # The total ground state energy is E_g = 4*omega + E_hop
    
    print("The ground state energy is given by the expression: E_g = 4*omega + E_hop")
    print("The hopping energy E_hop is calculated by summing the four lowest single-particle energies.")
    print("E_hop = -2*J * (cos(2*pi*m/7)) for the occupied states m = {0, 1, -1, 2}")
    print("\nSymbolically, the expression for E_hop is:")
    print("E_hop = -2 * J * (cos(0) + 2 * cos(2*pi/7) + cos(4*pi/7))\n")

    print("Let's calculate the numerical values for the terms in the equation:")
    print(f"The number of photons N_ph = {N_ph}")
    print(f"The constant energy term is {N_ph}*omega")
    
    print(f"\nThe numbers for the hopping term calculation are:")
    print(f"cos(0) = {cos0:.10f}")
    print(f"cos(2*pi/7) = {cos1:.10f}")
    print(f"cos(4*pi/7) = {cos2:.10f}")
    
    print(f"\nThe sum of cosines factor is: {cos0:.2f} + 2*({cos1:.10f}) + ({cos2:.10f}) = {sum_of_cosines:.10f}")
    
    print(f"\nSo, E_hop = -2 * J * ({sum_of_cosines:.10f})")
    print(f"E_hop = {j_coefficient:.10f} * J")

    print("\n-------------------------------------------------------------")
    print("The final expression for the ground state energy is:")
    print(f"E_g = {N_ph}*omega + ({j_coefficient:.10f})*J")
    print("-------------------------------------------------------------")

solve_hamiltonian()