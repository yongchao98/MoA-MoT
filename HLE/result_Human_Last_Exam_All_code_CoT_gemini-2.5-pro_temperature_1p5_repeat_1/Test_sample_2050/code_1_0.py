import numpy as np

def solve_cavity_problem():
    """
    Calculates the ground state energy for a system of 4 photons in 7 coupled
    cavities in a ring, in the limit of infinite on-site repulsion U.
    """
    # Parameters of the problem
    N = 7  # Number of cavities
    M = 4  # Number of photons

    # In the U -> infinity limit, the problem maps to M non-interacting fermions on an N-site ring.
    # The single-particle energies are E_m = -2*J*cos(2*pi*m/N).
    # We need to fill the M=4 lowest energy levels.
    # The quantum numbers 'm' for the 7 energy levels are m = 0, +/-1, +/-2, +/-3.
    # The lowest energy states correspond to m values that maximize cos(2*pi*m/N).
    # These are m=0, followed by m=+/-1, then m=+/-2, etc.
    # The M=4 lowest energy states correspond to m = 0, 1, -1, and one from m=+/-2.
    # We sum their energies to find the hopping contribution to the ground state energy.
    
    # E_hop = E(m=0) + E(m=1) + E(m=-1) + E(m=2)
    # E_hop/J = -2*cos(0) - 2*cos(2*pi/7) - 2*cos(-2*pi/7) - 2*cos(4*pi/7)
    # Since cos(-x) = cos(x), this simplifies to:
    # E_hop/J = -2 * (cos(0) + 2*cos(2*pi/7) + cos(4*pi/7))
    
    # The full ground state energy is E_ground = M*omega + E_hop

    print(f"The ground state energy E_ground for N={N} sites and M={M} photons is:")
    print(f"E_ground = {M}*omega + E_hopping")
    print("where E_hopping is the contribution from the hopping term.")
    print("\nThe hopping energy is the sum of the M lowest single-particle energies E_m = -2*J*cos(2*pi*m/N).")
    print(f"The occupied states are for m = 0, 1, -1, 2.")
    print("So, the hopping energy is: E_hopping = -2*J * (cos(2*pi*0/7) + 2*cos(2*pi*1/7) + cos(2*pi*2/7))")

    # Now we calculate the numerical values
    cos_0 = np.cos(0)
    cos_2pi_over_7 = np.cos(2 * np.pi / N)
    cos_4pi_over_7 = np.cos(4 * np.pi / N)

    # Calculate the sum inside the parenthesis
    coeff_sum = cos_0 + 2 * cos_2pi_over_7 + cos_4pi_over_7
    
    # Calculate the final coefficient for J
    j_coeff = -2 * coeff_sum

    print("\nLet's calculate the numerical value of the expression:")
    print(f"cos(0) is: {cos_0:.5f}")
    print(f"The number 2 appears because of the m=+1 and m=-1 degeneracy.")
    print(f"cos(2*pi/7) is: {cos_2pi_over_7:.5f}")
    print(f"cos(4*pi/7) is: {cos_4pi_over_7:.5f}")

    print("\nSubstituting these numbers into the equation for the hopping energy:")
    print(f"E_hopping = -2 * ({cos_0:.5f} + 2*{cos_2pi_over_7:.5f} + {cos_4pi_over_7:.5f}) * J")
    print(f"E_hopping = -2 * ({coeff_sum:.5f}) * J")
    print(f"E_hopping = ({j_coeff:.5f}) * J")
    
    print("\nTherefore, the final equation for the ground state energy is:")
    # Printing each number in the final equation
    print(f"{M} * omega + ({j_coeff:.5f}) * J")

solve_cavity_problem()