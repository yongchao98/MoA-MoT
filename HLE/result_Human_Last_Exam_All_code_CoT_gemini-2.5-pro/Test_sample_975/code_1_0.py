def display_magnetic_field_solution():
    """
    This function prints the derived expressions for the magnetic field H
    in the two specified regions, corresponding to the correct answer.
    """
    
    # Define symbolic representations for clarity in the print statements.
    M0 = "M_0"
    Rp = "R_p"
    R = "R"
    r = "r"
    cos_theta = "cos(θ)"
    sin_theta = "sin(θ)"
    ir = "î_r"
    itheta = "î_θ"

    print("Based on the derivation, the correct expressions for the magnetic field H are given by Option B.")
    print("The code below prints these final expressions in a clear format.\n")

    # --- Region 1: Inside the shield (0 < r < R_p) ---
    print("="*40)
    print(f"Region 1: 0 < r < {Rp}")
    print("="*40)
    
    # The full vector expression for H in Region 1
    # H = M_0 * (2*Rp^3 + R^3) / (3*R^3) * (-cos(θ)*î_r + sin(θ)*î_θ)
    
    coeff1_str = f"{M0} * (2*{Rp}^3 + {R}^3) / (3*{R}^3)"
    vector1_str = f"(-{cos_theta} * {ir} + {sin_theta} * {itheta})"
    
    print("The magnetic field H(r, θ) is:")
    print(f"H = ( {coeff1_str} ) * {vector1_str}\n")
    
    print("In component form:")
    # H_r component
    hr1_coeff_str = f"-({coeff1_str})"
    print(f"H_r = {hr1_coeff_str} * {cos_theta}")
    
    # H_theta component
    htheta1_coeff_str = f"({coeff1_str})"
    print(f"H_θ = {htheta1_coeff_str} * {sin_theta}")

    # --- Region 2: Between the shield and the conductor (R_p < r < R) ---
    print("\n" + "="*40)
    print(f"Region 2: {Rp} < r < {R}")
    print("="*40)

    # The components of H in Region 2
    # H_r = -2*M_0/3 * [ (R_p/R)^3 - (R_p/r)^3 ] * cos(θ)
    # H_θ = M_0/3 * [ 2*(R_p/R)^3 + (R_p/r)^3 ] * sin(θ)

    print("The magnetic field components H(r, θ) are:")
    
    # H_r component
    hr2_coeff_str = f"- (2*{M0}/3) * [({Rp}/{R})^3 - ({Rp}/{r})^3]"
    print(f"H_r = ( {hr2_coeff_str} ) * {cos_theta}")
    
    # H_theta component
    htheta2_coeff_str = f"( {M0}/3 ) * [2*({Rp}/{R})^3 + ({Rp}/{r})^3]"
    print(f"H_θ = ( {htheta2_coeff_str} ) * {sin_theta}")

if __name__ == '__main__':
    display_magnetic_field_solution()
