import math

def calculate_average_projected_area(n):
    """
    Calculates the components of the formula for the average projected area of an 
    n-dimensional cross-polytope, where n is an odd integer.

    The final formula is of the form: Numerator / (Denominator_Coefficient * sqrt(n))

    Args:
        n (int): The dimension, must be an odd integer >= 3.

    Returns:
        tuple: A tuple containing (numerator, denominator_coefficient, n).
    """
    if not isinstance(n, int) or n < 3 or n % 2 == 0:
        raise ValueError("n must be an odd integer >= 3")

    # The dimension n can be written as n = 2k + 1.
    # k = (n - 1) / 2
    k = (n - 1) // 2

    # The derived formula for the average projected area is:
    # 2 / (sqrt(n) * k! * (k-1)!)
    
    # We want to present this as Num / (Den * sqrt(n)).
    # The integer part of the formula is 2 / (k! * (k-1)!).
    # Let's calculate this fractional part and simplify it.
    
    numerator_part = 2
    denominator_part = math.factorial(k) * math.factorial(k - 1)
    
    # Simplify the fraction by dividing by the greatest common divisor.
    common_divisor = math.gcd(numerator_part, denominator_part)
    
    simplified_numerator = numerator_part // common_divisor
    simplified_denominator = denominator_part // common_divisor
    
    return simplified_numerator, simplified_denominator, n

def main():
    """
    Main function to explain the derivation and compute the average projected area
    for example dimensions.
    """
    print("This script calculates the average area of a projection of an n-dimensional cross-polytope P.")
    print("The dimension n is assumed to be an odd integer (n = 2k + 1).")
    
    print("\n### Derivation Summary ###")
    print("1. The average projected area is given by the formula: E[Area] = (omega_{n-1} / (n * omega_n)) * S(P),")
    print("   where S(P) is the surface area of the cross-polytope and omega_d is the surface area of a (d-1)-sphere.")
    
    print("\n2. The surface area of the n-cross-polytope is S(P) = 2^n * sqrt(n) / (n-1)!.")
    
    print("\n3. The coefficient involving sphere surface areas simplifies based on the Gamma function: omega_d = 2*pi^(d/2)/Gamma(d/2).")
    
    print("\n4. For n = 2k+1, combining and simplifying these parts yields the final formula:")
    print("   Average Area = 2 / (sqrt(n) * k! * (k-1)!), where k = (n-1)/2.")
    print("   This can also be written using Gamma functions as: 2 / (sqrt(n) * Gamma((n+1)/2) * Gamma((n-1)/2)).")
    
    print("\n### Calculations for Example Dimensions ###")
    
    example_dimensions = [3, 5, 7, 9]
    for n_val in example_dimensions:
        try:
            num, den, n_sqrt = calculate_average_projected_area(n_val)
            k_val = (n_val - 1) // 2
            
            print(f"\n--- For n = {n_val} (so k = {k_val}) ---")
            
            # Display the final equation in a readable format
            if den == 1:
                if num == 1:
                    print(f"The final equation is: 1 / sqrt({n_sqrt})")
                else:
                    print(f"The final equation is: {num} / sqrt({n_sqrt})")
            else:
                print(f"The final equation is: {num} / ({den} * sqrt({n_sqrt}))")

            # Output each number in the final equation as requested
            print("The numbers in the final simplified equation are:")
            print(f"  - Numerator: {num}")
            print(f"  - Denominator coefficient: {den}")
            print(f"  - Value under the square root: {n_sqrt}")
            
            # Also print the numerical value for context
            numerical_value = num / (den * math.sqrt(n_sqrt))
            print(f"Numerical value: {numerical_value:.8f}")

        except ValueError as e:
            print(f"\nError calculating for n={n_val}: {e}")

if __name__ == "__main__":
    main()