def display_magnetic_field_solution():
    """
    This function presents the solution to the magnetostatics problem
    by printing the final, derived equations for the magnetic field H.
    The equations correspond to the correct answer choice.
    The numerical factors in the equations are explicitly printed to
    satisfy the problem's output requirements.
    """

    print("Based on solving the magnetostatic boundary value problem, the derived magnetic field H(r, theta) is presented below.")
    print("The solution corresponds to choice B.")

    print("\n" + "="*60 + "\n")
    print("--- In the region 0 < r < Rp (within the shield) ---")
    print("The magnetic field is uniform and given by:")
    print("H = C * (-cos(theta) * i_r + sin(theta) * i_theta)")
    print("where the coefficient C is:")
    print("C = M0 * (2*Rp^3 + R^3) / (3*R^3)")
    print("\nTo satisfy the instruction 'output each number in the final equation':")
    print("Equation form: H = M0 * ( A*Rp^3 + B*R^3 ) / ( C*R^3 ) * ( D*cos(theta)*i_r + E*sin(theta)*i_theta )")
    print(f"The numbers are: A = 2, B = 1, C = 3, D = -1, E = 1")

    print("\n" + "="*60 + "\n")

    print("--- In the region Rp < r < R (between shield and conductor) ---")
    print("The magnetic field is non-uniform and has two components (radial and tangential):")
    
    print("\nRadial component (H_r):")
    print("H_r = Coeff_r * M0 * [ (Rp/R)^3 - (Rp/r)^3 ] * cos(theta)")
    print("The coefficient Coeff_r is -2/3.")
    
    print("\nTangential component (H_theta):")
    print("H_theta = Coeff_theta * M0 * [ 2*(Rp/R)^3 + (Rp/r)^3 ] * sin(theta)")
    print("The coefficient Coeff_theta is 1/3.")

    print("\nTo satisfy the instruction 'output each number in the final equation':")
    print("Equation form for H_r: H_r = (A/B) * M0 * [ (Rp/R)^3 - (Rp/r)^3 ] * cos(theta)")
    print(f"The numbers are: A = -2, B = 3")
    print("\nEquation form for H_theta: H_theta = (C/D) * M0 * [ E*(Rp/R)^3 + F*(Rp/r)^3 ] * sin(theta)")
    print(f"The numbers are: C = 1, D = 3, E = 2, F = 1")


if __name__ == '__main__':
    display_magnetic_field_solution()
