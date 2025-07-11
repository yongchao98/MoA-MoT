import math

def evaluate_1s_slater_integral():
    """
    Calculates the integral <phi_i| 1/r |phi_j> for 1s Slater orbitals.
    """
    try:
        # Get user input for the orbital exponents
        zeta_i_str = input("Enter the value for orbital exponent zeta_i: ")
        zeta_i = float(zeta_i_str)
        
        zeta_j_str = input("Enter the value for orbital exponent zeta_j: ")
        zeta_j = float(zeta_j_str)

        if zeta_i <= 0 or zeta_j <= 0:
            print("Error: Orbital exponents must be positive numbers.")
            return

    except ValueError:
        print("Invalid input. Please enter valid numbers.")
        return

    print("-" * 50)
    print(f"Calculating the integral for zeta_i = {zeta_i} and zeta_j = {zeta_j}")
    print("-" * 50)
    
    # Check for the special case where i = j
    if zeta_i == zeta_j:
        zeta = zeta_i
        result = zeta
        print("This is the special case where zeta_i = zeta_j = zeta.")
        print("The integral <phi|1/r|phi> simplifies to the expectation value of 1/r.")
        print("The final equation is: <phi|1/r|phi> = zeta")
        print(f"Result = {result}")

    else:
        # General case for i != j
        # The analytical formula is I = (4 * (zeta_i * zeta_j)^(3/2)) / (zeta_i + zeta_j)^2
        print("The general formula for the integral <phi_i|1/r|phi_j> is:")
        print("I = (4 * (zeta_i * zeta_j)**(3/2)) / (zeta_i + zeta_j)**2\n")

        print("Plugging in the values:")
        # Step-by-step calculation output
        term1_val = zeta_i * zeta_j
        print(f"I = (4 * ({zeta_i} * {zeta_j})**1.5) / ({zeta_i} + {zeta_j})**2")
        
        term1_pow_val = term1_val ** 1.5
        print(f"I = (4 * {term1_val}**1.5) / ({zeta_i + zeta_j})**2")
        
        numerator = 4 * term1_pow_val
        print(f"I = (4 * {term1_pow_val}) / ({zeta_i + zeta_j})**2")
        
        denominator = (zeta_i + zeta_j) ** 2
        print(f"I = {numerator} / {denominator}")
        
        result = numerator / denominator
        print(f"\nFinal Result: I = {result}")

if __name__ == '__main__':
    evaluate_1s_slater_integral()