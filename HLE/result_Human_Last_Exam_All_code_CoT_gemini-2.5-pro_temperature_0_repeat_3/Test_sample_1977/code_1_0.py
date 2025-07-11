import math

def solve_norm():
    """
    Calculates the 1-norm of the correlation matrix T for the state J_n for a given even n.
    """
    try:
        n_str = input("Enter a non-negative even integer n: ")
        n = int(n_str)
        if n < 0 or n % 2 != 0:
            print("Error: n must be a non-negative even integer.")
            return
    except ValueError:
        print("Error: Invalid input. Please enter an integer.")
        return

    # The derived formula for the 1-norm of the correlation matrix T is:
    # ||T||_1 = (2**(n+1) / (1 + 3**n)) * (3**(n+1) + 2 * 6**n - 2**(n+1) - 1)
    # This formula is valid for any even integer n >= 0.

    # Calculate the terms in the formula
    
    # Denominator of the prefactor
    denominator = 1 + 3**n
    
    # Numerator of the prefactor
    prefactor_num = 2**(n + 1)

    # Terms inside the parenthesis (Sigma)
    term1_sigma = 3**(n + 1)
    term2_sigma = 2 * (6**n)
    term3_sigma = 2**(n + 1)
    term4_sigma = 1
    
    sigma = term1_sigma + term2_sigma - term3_sigma - term4_sigma
    
    # Full numerator
    numerator = prefactor_num * sigma
    
    # Final result
    # Use integer division if possible, otherwise float division.
    if numerator % denominator == 0:
        result = numerator // denominator
    else:
        result = numerator / denominator

    print("\n" + "="*40)
    print(f"Calculation for n = {n}")
    print("="*40)
    
    print("\nThe 1-norm of the correlation matrix T is given by the formula:")
    print(f"||T||_1 = (2^(n+1) / (1 + 3^n)) * (3^(n+1) + 2*6^n - 2^(n+1) - 1)")
    
    print("\nStep 1: Calculate the main components of the formula.")
    
    print(f"  - Prefactor: 2^({n}+1) / (1 + 3^{n}) = {prefactor_num} / {denominator}")
    
    print(f"  - Sum term (Sigma): 3^({n}+1) + 2*6^{n} - 2^({n}+1) - 1")
    print(f"    = {term1_sigma} + {term2_sigma} - {term3_sigma} - {term4_sigma}")
    print(f"    = {sigma}")

    print("\nStep 2: Combine the components to get the final equation.")
    print("The final equation with the calculated numbers is:")
    print(f"||T||_1 = ({prefactor_num} * {sigma}) / {denominator}")
    print(f"        = {numerator} / {denominator}")
    print(f"        = {result}")
    print("="*40)


if __name__ == '__main__':
    solve_norm()
