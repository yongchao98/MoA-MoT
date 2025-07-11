import numpy as np

def calculate_ground_state_energy():
    """
    Calculates the ground state energy for the given Hamiltonian.
    """
    N = 7
    num_photons = 4

    # The ground state energy for the U->inf Hubbard-like model
    # is found by summing the lowest single-particle energies of the hopping term
    # on a ring, plus a constant offset from the omega term.
    
    # Single-particle energies are E_m = -2J * cos(2*pi*m/N).
    # We need to sum the energies for the 4 modes with the lowest energy.
    # The lowest energies correspond to the largest values of cos(2*pi*m/N).
    
    # For N=7, the unique cosine values in decreasing order are for m=0, 1, 2, 3.
    # cos(0) - non-degenerate
    # cos(2*pi/7) = cos(12*pi/7) (for m=1,6) - degenerate
    # cos(4*pi/7) = cos(10*pi/7) (for m=2,5) - degenerate
    
    # The 4 lowest energy states correspond to filling modes m=0, m=1, m=6, and m=2.
    
    # The hopping energy is E_hop = -2J * (cos(0) + cos(2*pi/7) + cos(12*pi/7) + cos(4*pi/7))
    # E_hop = -2J * (cos(0) + 2*cos(2*pi/7) + cos(4*pi/7))
    
    print("Step 1: The problem simplifies to 4 non-interacting hard-core bosons (fermions) on a 7-site ring.")
    print("The ground state energy is E_gs = 4\u03C9 + E_hop.")

    print("\nStep 2: Calculate the hopping part of the energy, E_hop.")
    print("E_hop is the sum of the 4 lowest single-particle energies.")
    print("E_hop = -2*J * (cos(0) + 2*cos(2\u03C0/7) + cos(4\u03C0/7))")
    
    # Numerically calculate the terms in the coefficient
    c0 = np.cos(0)
    c1 = np.cos(2 * np.pi / N)
    c2 = np.cos(4 * np.pi / N)
    
    total_coeff_sum = c0 + 2*c1 + c2
    C_J = 2 * total_coeff_sum

    print("\nStep 3: Substitute the numerical values into the equation.")
    print("The final equation for the ground state energy is E_gs = 4\u03C9 - C_J * J, where C_J is a numerical coefficient.")
    print(f"C_J = 2 * (cos(0) + 2*cos(2\u03C0/7) + cos(4\u03C0/7))")
    print(f"C_J = 2 * ({c0:.5f} + 2*{c1:.5f} + {c2:.5f})")
    print(f"C_J = 2 * ({total_coeff_sum:.5f})")
    print(f"C_J = {C_J:.5f}")

    print(f"\nFinal Result: The ground state energy is E_gs = {num_photons}\u03C9 - {C_J}*J")

    # This is the numerical answer expected for the coefficient C_J.
    return C_J

if __name__ == '__main__':
    C_J = calculate_ground_state_energy()
    # The problem asks to calculate the ground state energy. Since omega and J are not given,
    # the answer is an expression. The non-trivial calculation is the numerical coefficient.
    # So we provide the coefficient C_J = 2*(1 + 2*cos(2*pi/7) + cos(4*pi/7)).
    print(f"\nThe ground state energy is E_gs = 4*omega - {C_J}*J.")
    print("The numerical coefficient of (-J) is presented as the final answer.")

# The final answer is the coefficient C_J.
final_answer = 2 * (np.cos(0) + 2 * np.cos(2 * np.pi / 7) + np.cos(4 * np.pi / 7))
# print(f"<<<{final_answer}>>>") # Suppressing this print as it shouldn't be in the code block itself