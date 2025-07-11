import math

def solve_energy_minimization():
    """
    This function calculates the minimum total heat energy based on the derived analytical solution.
    The steps are outlined in the explanation.
    """

    print("Step 1: The problem reduces to finding the minimum of E = 1 + I_1/2,")
    print("subject to the constraint I_1^2 - (4/3)*I_1 - 32 <= 0, where I_1 is an integral.")
    print("\nStep 2: We solve the quadratic equation 3*y^2 - 4*y - 96 = 0 to find the bounds for I_1.")

    # Coefficients of the quadratic equation 3*y^2 - 4*y - 96 = 0
    a = 3
    b = -4
    c = -96

    # Calculate the discriminant
    discriminant = b**2 - 4*a*c
    
    print(f"\nThe coefficients are a={a}, b={b}, c={c}.")
    print(f"The discriminant is b^2 - 4ac = ({b})^2 - 4*({a})*({c}) = {discriminant}.")

    # Calculate the two roots
    sqrt_discriminant = math.sqrt(discriminant)
    root1 = ( -b + sqrt_discriminant ) / (2*a)
    root2 = ( -b - sqrt_discriminant ) / (2*a)
    
    print(f"The roots for I_1 are ({ -b } + sqrt({discriminant})) / (2*{a}) and ({ -b } - sqrt({discriminant})) / (2*{a}).")

    # The minimum value for I_1 is the smaller of the two roots.
    I1_min = min(root1, root2)
    print(f"\nStep 3: The allowed range for I_1 is [{min(root1, root2):.4f}, {max(root1, root2):.4f}].")
    print(f"To minimize E, we take the minimum possible value for I_1, which is {I1_min:.4f}.")
    
    # Calculate the minimum energy
    E_min = 1 + I1_min / 2
    
    print("\nStep 4: The minimum energy is E_min = 1 + I_1_min / 2.")
    print("The final equation for the minimum energy is E_min = (4 - sqrt(73)) / 3.")
    
    # Let's show the symbolic calculation
    # discriminant = 1168 = 16 * 73
    # sqrt(1168) = 4 * sqrt(73)
    # I1_min = (4 - 4*sqrt(73))/6 = (2 - 2*sqrt(73))/3
    # E_min = 1 + (1/2) * (2 - 2*sqrt(73))/3 = 1 + (1-sqrt(73))/3 = (3+1-sqrt(73))/3 = (4-sqrt(73))/3
    num_4 = 4
    num_73 = 73
    num_3 = 3
    print(f"Final symbolic answer: E_min = ({num_4} - sqrt({num_73})) / {num_3}")

    print(f"\nNumerical value of the minimum energy E_min is: {E_min}")

solve_energy_minimization()