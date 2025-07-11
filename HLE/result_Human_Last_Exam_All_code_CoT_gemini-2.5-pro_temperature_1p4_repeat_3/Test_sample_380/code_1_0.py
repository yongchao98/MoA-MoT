import numpy as np

def calculate_helicity_comparison():
    """
    Demonstrates the difference between calculating a helicity-related term
    in a simple Radial-Tangential-Normal (RTN) system versus a rigorous
    Field-Aligned Coordinate (FAC) system.

    We use a proxy for helicity: Im(b_perp1^* * b_perp2), where b_perp
    are the complex magnetic field components of the wave perpendicular
    to the main field direction.
    """

    # 1. Define the background magnetic field (B0) at L1 in RTN coordinates (nT).
    # R (Radial), T (Tangential), N (Normal).
    # A typical Parker spiral at 1 AU has Br and Bt components of similar magnitude.
    # We choose Br = -4 nT (points towards Sun) and Bt = 4 nT. Bn is small.
    B0_rtn = np.array([-4.0, 4.0, 0.5])
    B0_magnitude = np.linalg.norm(B0_rtn)
    
    # 2. Define a sample wave magnetic field (delta_b) in RTN coordinates (nT).
    # The components are complex numbers to represent wave phase.
    # For AIC waves, fluctuations are primarily perpendicular to B0.
    delta_b_rtn = np.array([0.1 + 0.05j, 0.5 + 0.2j, 0.3 - 0.4j])
    
    print("--- Setup ---")
    print(f"Background Magnetic Field B0 (RTN): [{B0_rtn[0]:.2f}, {B0_rtn[1]:.2f}, {B0_rtn[2]:.2f}] nT")
    print(f"Wave Fluctuation delta_b (RTN): [({delta_b_rtn[0].real:.2f}{delta_b_rtn[0].imag:+.2f}j), "
          f"({delta_b_rtn[1].real:.2f}{delta_b_rtn[1].imag:+.2f}j), "
          f"({delta_b_rtn[2].real:.2f}{delta_b_rtn[2].imag:+.2f}j)] nT")
    print("-" * 20 + "\n")

    # 3. Approximate Calculation (assuming B0 is radial)
    # The perpendicular components are assumed to be T and N.
    # We calculate Im(b_T* * b_N)
    b_T = delta_b_rtn[1]
    b_N = delta_b_rtn[2]
    helicity_term_approx = np.imag(np.conj(b_T) * b_N)

    print("--- 1. Approximate Calculation (Assuming B0 is Radial) ---")
    print("The perpendicular components are assumed to be the Tangential (T) and Normal (N) components.")
    # We use f-string formatting to explicitly show the numbers in the equation.
    print(f"Equation: Im(b_T* * b_N)")
    print(f"  b_T = ({b_T.real:.2f}{b_T.imag:+.2f}j), b_N = ({b_N.real:.2f}{b_N.imag:+.2f}j)")
    print(f"Result: {helicity_term_approx:.4f}\n")
    
    # 4. Rigorous Calculation (using a Field-Aligned Coordinate system)
    # Define the FAC basis vectors from B0.
    z_fac = B0_rtn / B0_magnitude  # Parallel to B0
    # Create a y-axis perpendicular to the B0-R plane. Let R_hat be [1, 0, 0]
    y_fac = np.cross(z_fac, np.array([1.0, 0.0, 0.0]))
    y_fac = y_fac / np.linalg.norm(y_fac) # Normalize
    # The x-axis completes the right-handed system.
    x_fac = np.cross(y_fac, z_fac)

    # Create the transformation matrix from FAC to RTN
    # The columns are the new basis vectors expressed in the old system.
    transform_matrix = np.array([x_fac, y_fac, z_fac]).T
    # To go from RTN to FAC, we need the inverse, which is the transpose for an orthonormal matrix.
    inverse_transform = transform_matrix.T
    
    # Transform the wave vector into the FAC system.
    delta_b_fac = inverse_transform @ delta_b_rtn

    # In FAC, the perpendicular components are the first two (x_fac and y_fac).
    b_perp1 = delta_b_fac[0]
    b_perp2 = delta_b_fac[1]
    helicity_term_rigorous = np.imag(np.conj(b_perp1) * b_perp2)
    
    print("--- 2. Rigorous Calculation (Using Field-Aligned Coordinates) ---")
    print("The perpendicular components are properly calculated in a coordinate system aligned with B0.")
    print(f"Equation: Im(b_perp1* * b_perp2)")
    print(f"  b_perp1 = ({b_perp1.real:.2f}{b_perp1.imag:+.2f}j), b_perp2 = ({b_perp2.real:.2f}{b_perp2.imag:+.2f}j)")
    print(f"Result: {helicity_term_rigorous:.4f}\n")
    
    # 5. Compare the results
    difference = abs(helicity_term_approx - helicity_term_rigorous)
    print("--- Comparison ---")
    print(f"Approximate Value: {helicity_term_approx:.4f}")
    print(f"Rigorous Value:    {helicity_term_rigorous:.4f}")
    print(f"Absolute Difference: {difference:.4f}")

if __name__ == '__main__':
    calculate_helicity_comparison()