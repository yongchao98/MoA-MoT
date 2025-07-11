import math

def calculate_ground_state_energy():
    """
    This function calculates the ground state energy for a system of 4 photons
    in 7 coupled nonlinear optical cavities in the hard-core limit (U -> infinity).
    
    The plan is as follows:
    1.  The system is modeled as 4 non-interacting fermions on a 7-site ring due to
        the U -> infinity limit (hard-core bosons -> fermions).
    2.  The single-particle energies are epsilon_m = omega - 2*J*cos(2*pi*m/7).
    3.  The ground state energy is the sum of the 4 lowest single-particle energies.
    4.  These correspond to filling the momentum states m=0, m=1, m=6, and m=2.
    5.  The total energy is E_GS = 4*omega - 2*J*(cos(0) + cos(2pi/7) + cos(12pi/7) + cos(4pi/7)).
    6.  This simplifies to E_GS = 4*omega - 2*J*(1 + 2*cos(2pi/7) + cos(4pi/7)).
    7.  This script calculates the numerical coefficient for J and prints the final expression.
    """
    
    # System parameters
    N = 7  # Number of cavities
    P = 4  # Number of photons

    # The coefficient of omega is simply the number of photons.
    omega_coeff = float(P)

    # --- Calculation for the coefficient of J ---
    
    # The term comes from -2 * sum of cosines for the 4 lowest energy modes.
    # The sum is (cos(0) + cos(2pi/7) + cos(12pi/7) + cos(4pi/7)).
    # Since cos(12pi/7) = cos(2pi/7), this is (1 + 2*cos(2pi/7) + cos(4pi/7)).

    # Let's define the components for the calculation.
    constant_term = 1.0  # from cos(0)
    cos_factor = 2.0     # factor for the cos(2pi/7) term

    # Calculate the required cosine values
    angle1 = 2 * math.pi / N
    cos_val1 = math.cos(angle1)

    angle2 = 4 * math.pi / N
    cos_val2 = math.cos(angle2)

    # Calculate the final coefficient for J
    J_coeff = -2.0 * (constant_term + cos_factor * cos_val1 + cos_val2)

    # --- Output the results ---
    
    print("The ground state energy E_GS can be expressed in the form: E_GS = A * omega + B * J")
    print("\nHere are the calculated values for the coefficients A and B:")
    
    # Output the value for A
    print(f"\nCoefficient A (the factor of omega) = {omega_coeff:.1f}")

    # Output the components and final value for B
    print("\nCoefficient B (the factor of J) is derived from: B = -2 * (1 + 2*cos(2*pi/7) + cos(4*pi/7))")
    print("The numerical values used are:")
    print(f"  cos(2*pi/7) = {cos_val1:.6f}")
    print(f"  cos(4*pi/7) = {cos_val2:.6f}")
    print(f"Substituting these values gives B = -2 * ({constant_term:.1f} + {cos_factor:.1f} * {cos_val1:.6f} + {cos_val2:.6f})")
    print(f"Resulting in B = {J_coeff:.6f}")

    print("\n------------------------------------------------------------------")
    print("The final expression for the ground state energy is therefore:")
    print(f"E_GS = {omega_coeff:.1f} * omega + ({J_coeff:.6f}) * J")
    print("------------------------------------------------------------------")


calculate_ground_state_energy()
