import math

def print_final_answer():
    """
    This function prints the final, symbolic equations for the magnetic field H
    in the two specified regions, as determined by the derivation.
    """
    print("Based on the derivation, the correct expressions are given in Option B.")
    print("\nThe final equations for the magnetic field H(r, theta) are:")

    print("\n--- In the region 0 < r < R_p ---")
    # Equation for H inside the shield
    # H = M0 * (2*Rp^3 + R^3) / (3*R^3) * (-cos(theta)*i_r + sin(theta)*i_theta)
    print("The vector equation for H is:")
    print("H = (M_0 * ((2 * (R_p**3)) + R**3)) / (3 * (R**3)) * ((-1 * cos(theta) * i_r) + (1 * sin(theta) * i_theta))")

    print("\n--- In the region R_p < r < R ---")
    # Equations for H components between the shield and the conductor
    # H_r = -(2*M0/3) * [ (Rp/R)^3 - (Rp/r)^3 ] * cos(theta)
    # H_theta = (M0/3) * [ 2*(Rp/R)^3 + (Rp/r)^3 ] * sin(theta)
    print("The vector components of H are:")
    print("H_r = (-(2 * M_0 / 3) * (((R_p/R)**3) - ((R_p/r)**3))) * cos(theta)")
    print("H_theta = ((1 * M_0 / 3) * ((2 * ((R_p/R)**3)) + ((R_p/r)**3))) * sin(theta)")

# Execute the function to print the result
print_final_answer()