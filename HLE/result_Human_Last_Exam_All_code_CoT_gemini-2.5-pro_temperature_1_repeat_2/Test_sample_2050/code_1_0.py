import math

def calculate_ground_state_energy():
    """
    Calculates the ground state energy for the given Hamiltonian.

    The Hamiltonian is:
    H = Sum_i (omega * a_i^dag a_i + U/2 * a_i^dag a_i^dag a_i a_i - J * (a_i^dag a_{i+1} + a_i^dag a_{i-1}))
    with N=7 cavities in a ring and 4 photons.

    In the limit U -> infinity, no two photons can occupy the same site.
    This is the hard-core boson limit, which is equivalent to non-interacting spinless fermions in 1D.
    The problem reduces to filling the lowest single-particle energy levels of a tight-binding ring.

    The single-particle energies are E_m = -2J * cos(2*pi*m/N), with m = 0, +/-1, +/-2, ...
    For N=7, the lowest energies correspond to m values that maximize cos(2*pi*m/7).
    The ordering of energy levels from lowest to highest is for m = 0, {1, -1}, {2, -2}, {3, -3}.

    For 4 photons (fermions), we fill the 4 lowest energy states:
    - One particle in the m=0 state.
    - Two particles in the degenerate m=+1 and m=-1 states.
    - One particle in one of the degenerate m=+2 or m=-2 states.

    The total energy is the sum of these single-particle energies plus the constant term.
    E_ground = sum(E_single_particle) + constant_term
    The constant term from on-site energy is Sum_i(omega * n_i) = 4 * omega.
    The hopping energy is Sum_{filled states} E_m.
    E_hopping = E_0 + E_1 + E_{-1} + E_2
               = -2J*(cos(0) + cos(2*pi/7) + cos(-2*pi/7) + cos(4*pi/7))
               = -2J*(1 + 2*cos(2*pi/7) + cos(4*pi/7))
    """
    # System parameters
    N = 7  # Number of cavities
    n_photons = 4  # Number of photons

    # The coefficient of omega is the total number of photons
    C_omega = n_photons

    # The coefficient of J is determined by the sum of cosines from the occupied energy levels
    # The term inside the parenthesis is sum_{occupied k} 2*cos(k)
    # Here we calculate sum_{occupied m} cos(2*pi*m/N)
    # The occupied states are m=0, 1, -1, 2
    sum_of_cosines = (math.cos(2 * math.pi * 0 / N) +
                      math.cos(2 * math.pi * 1 / N) +
                      math.cos(2 * math.pi * -1 / N) +
                      math.cos(2 * math.pi * 2 / N))

    C_J = -2 * sum_of_cosines

    # Print the final equation with all the numbers
    print("The ground state energy E_ground is given by the equation:")
    print(f"E_ground = {C_omega} * omega + {C_J} * J")
    
    # As requested, printing each numerical coefficient in the final equation
    print("\nBreaking down the numbers in the final equation:")
    print(f"The number of photons is {n_photons}, which is the coefficient of omega.")
    print(f"The calculated coefficient of J is: {C_J}")

calculate_ground_state_energy()