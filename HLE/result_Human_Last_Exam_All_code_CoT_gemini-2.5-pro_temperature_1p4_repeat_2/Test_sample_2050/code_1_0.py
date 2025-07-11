import math

def calculate_ground_state_energy():
    """
    Calculates the ground state energy for the given system.

    The problem involves 4 photons in 7 coupled cavities in a ring,
    in the limit of infinite on-site repulsion (U -> infinity).
    This corresponds to a system of 4 hard-core bosons, which in 1D
    can be mapped to non-interacting spinless fermions.

    The ground state energy is the sum of the on-site energy part
    and the hopping energy part.
    - On-site energy: 4*omega, since there are 4 photons.
    - Hopping energy: Sum of the 4 lowest single-particle energies of
      a particle on a 7-site ring.
    
    The single-particle energies are given by epsilon_m = -2*J*cos(2*pi*m/N),
    where N=7 and m = 0, 1, ..., 6.
    
    We fill the energy levels from the bottom up:
    - 1 particle in the m=0 state.
    - 2 particles in the degenerate m=1, m=6 states.
    - 1 particle in one of the degenerate m=2, m=5 states.
    
    The code calculates the coefficient of J in the final energy expression.
    """
    N = 7
    P = 4

    # Sum of the cosine terms for the 4 lowest energy states (m=0, 1, 6, 2)
    # The term for J is -2 * J * sum(cos(k_m))
    cos_sum = math.cos(0) + math.cos(2 * math.pi / N) + math.cos(2 * math.pi * 6 / N) + math.cos(2 * math.pi * 2 / N)
    
    # cos(2*pi*6/N) is equal to cos(2*pi/N), so the sum is 1 + 2*cos(2*pi/7) + cos(4*pi/7)
    
    # The coefficient of J in the energy expression
    j_coefficient = -2 * cos_sum

    print("The ground state energy E_ground is given by the expression:")
    print(f"E_ground = {P}*omega + Hopping_Energy")
    print("\nThe symbolic form of the equation is:")
    print(f"E_ground = {P}*omega - 2*J*(1 + 2*cos(2*pi/{N}) + cos(4*pi/{N}))")
    
    print("\nNumerically, the equation is approximately:")
    print(f"E_ground = {P}*omega + ({j_coefficient:.10f})*J")

if __name__ == '__main__':
    calculate_ground_state_energy()
    # The value for the final answer tag
    N = 7
    cos_sum = 1 + 2 * math.cos(2 * math.pi / N) + math.cos(4 * math.pi / N)
    j_coefficient = -2 * cos_sum
    # print(f"\n<<<{j_coefficient:.10f}>>>")