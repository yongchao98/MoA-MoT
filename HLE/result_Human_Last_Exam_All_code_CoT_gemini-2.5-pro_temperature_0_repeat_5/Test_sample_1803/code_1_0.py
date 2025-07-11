import sys

# This script calculates the ground state energy shift for two interacting
# quantum harmonic oscillators using perturbation theory.

# Suppress the creation of __pycache__ directory
sys.dont_write_bytecode = True

def calculate_energy_shift():
    """
    This function prints the step-by-step derivation of the ground state
    energy shift (London dispersion force) between two quantum harmonic oscillators.
    """
    print("Step-by-step derivation of the ground state energy shift:")
    print("=========================================================")

    print("\nStep 1: The Interaction Hamiltonian (Perturbation H')")
    print("The system consists of two 3D quantum harmonic oscillators (QHOs), separated by a large distance R.")
    print("Each QHO consists of a positive and negative charge (+e, -e), forming an electric dipole.")
    print("The dipole moment of oscillator 'i' is p_i = e * x_i, where x_i is the displacement vector.")
    print("We align the system such that the separation vector R is along the z-axis.")
    print("For large R, the leading term in the Coulomb interaction is the dipole-dipole potential, which is our perturbation H':")
    print("H' = (1 / (4*pi*epsilon_0)) * [ (p1 . p2) / R^3 - 3 * (p1 . R_hat) * (p2 . R_hat) / R^3 ]")
    print("Substituting p_i = e*x_i and R_hat = z_hat, this becomes:")
    print("H' = (e^2 / (4*pi*epsilon_0 * R^3)) * [x1 . x2 - 3*(x1_z)*(x2_z)]")
    print("H' = (e^2 / (4*pi*epsilon_0 * R^3)) * [x1_x*x2_x + x1_y*x2_y - 2*x1_z*x2_z]")

    print("\nStep 2: Applying Perturbation Theory")
    print("The unperturbed ground state is |0> = |0>_1 * |0>_2, where |0>_i is the ground state of the i-th 3D QHO.")
    print("The first-order energy correction Delta_E^(1) = <0|H'|0> is zero.")
    print("This is because H' is linear in position operators x_i, and the ground state expectation value <0|x_i|0> is zero.")

    print("\nStep 3: Second-Order Energy Correction")
    print("The leading non-zero contribution is the second-order correction:")
    print("Delta_E = sum_{n!=0} |<n|H'|0>|^2 / (E0 - En)")
    print("\n- The Matrix Elements <n|H'|0>:")
    print("  The position operator x can be written as x = sqrt(hbar / (2*m*omega_0)) * (a + a_dagger).")
    print("  H' connects the ground state |0> only to states |n> where both oscillators are excited from the ground state to the first excited state.")
    print("  The non-zero matrix elements are for states like |1_x>_1|1_x>_2, |1_y>_1|1_y>_2, and |1_z>_1|1_z>_2.")
    print("  Let C = e^2 / (4*pi*epsilon_0 * R^3). The matrix elements are:")
    print("  <1x,1x|H'|0> = C * (hbar/(2*m*omega_0))")
    print("  <1y,1y|H'|0> = C * (hbar/(2*m*omega_0))")
    print("  <1z,1z|H'|0> = -2 * C * (hbar/(2*m*omega_0))")
    print("\n- The Energy Denominator (E0 - En):")
    print("  The ground state energy is E0 = (3/2)*hbar*omega_0 + (3/2)*hbar*omega_0 = 3*hbar*omega_0.")
    print("  The energy of an excited state like |1x,1x> is En = (5/2)*hbar*omega_0 + (5/2)*hbar*omega_0 = 5*hbar*omega_0.")
    print("  The energy denominator is E0 - En = 3*hbar*omega_0 - 5*hbar*omega_0 = -2*hbar*omega_0.")

    print("\nStep 4: Final Calculation")
    print("Summing the contributions from the three excited states:")
    print("Delta_E = (1 / (-2*hbar*omega_0)) * [ |<1x,1x|H'|0>|^2 + |<1y,1y|H'|0>|^2 + |<1z,1z|H'|0>|^2 ]")
    print("Delta_E = (1 / (-2*hbar*omega_0)) * [C^2 * (hbar/(2*m*omega_0))^2 * (1^2 + 1^2 + (-2)^2)]")
    print("Delta_E = (1 / (-2*hbar*omega_0)) * [C^2 * (hbar/(2*m*omega_0))^2 * 6]")
    print("Delta_E = -3 * C^2 * (hbar/(2*m*omega_0))^2 / (hbar*omega_0)")
    print("Substituting C = e^2 / (4*pi*epsilon_0 * R^3) and simplifying:")
    print("Delta_E = -3 * (e^2 / (4*pi*epsilon_0 * R^3))^2 * hbar^2 / (4*m^2*omega_0^2 * hbar*omega_0)")

    print("\n----------------------------------------------------")
    print("Final Result: The Ground State Energy Shift")
    print("----------------------------------------------------")
    final_expression = "- (3/4) * (e^4 * hbar) / ((4*pi*epsilon_0)^2 * m^2 * omega_0^3 * R^6)"
    print(f"Delta_E = {final_expression}")

    print("\nBreaking down the final equation as requested:")
    numerator_factor = 3
    denominator_factor = 4
    print(f"The numerator of the numerical factor is: {numerator_factor}")
    print(f"The denominator of the numerical factor is: {denominator_factor}")
    print("The other factors in the equation are:")
    print("  - In the numerator: e^4 (electron charge to the fourth power), hbar (reduced Planck constant)")
    print("  - In the denominator: (4*pi*epsilon_0)^2 (Coulomb constant factor squared), m^2 (mass squared), omega_0^3 (frequency cubed), R^6 (distance to the sixth power)")

if __name__ == '__main__':
    calculate_energy_shift()