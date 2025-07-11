def calculate_norm_for_Jn():
    """
    Calculates the 1-norm of the correlation matrix T for the state J_n for a given even n.
    """
    print("This script calculates the 1-norm of the correlation matrix for the quantum state J_n.")
    print("The formula is valid for non-negative even integers n.")
    
    try:
        # You can change the value of n here to any non-negative even integer.
        n_input = input("Please enter a non-negative even integer for n (e.g., 0, 2, 4, ...): ")
        n = int(n_input)
        if n < 0 or n % 2 != 0:
            print(f"Error: The provided number n = {n} is not a non-negative even integer.")
            return
    except (ValueError, TypeError):
        print("Error: Invalid input. Please enter a valid integer.")
        return

    # Numerator calculation using Python's arbitrary-precision integers
    # Numerator = 2**(n+1) * 3**n - 2**(n+1) + 3**(n+1) - 1
    term1_num = 2**(n + 1) * 3**n
    term2_num = 2**(n + 1)
    term3_num = 3**(n + 1)
    term4_num = 1
    numerator = term1_num - term2_num + term3_num - term4_num

    # Denominator calculation
    # Denominator = 2**(n+1) * (1 + 3**n)
    term1_den = 2**(n + 1)
    term2_den = 1 + 3**n
    denominator = term1_den * term2_den

    # Output the equation with the specific numbers for the given n
    print("\n" + "="*50)
    print(f"For n = {n}, the calculation is as follows:")
    print("="*50)
    print("Formula: ||T||_1 = (2^(n+1)*3^n - 2^(n+1) + 3^(n+1) - 1) / (2^(n+1)*(1 + 3^n))")
    print("\nSubstituting n with the given value:")
    
    # Printing the numerator calculation
    print(f"Numerator   = (2^({n+1}) * 3^{n}) - 2^({n+1}) + 3^({n+1}) - 1")
    print(f"            = ({term1_num}) - {term2_num} + {term3_num} - {term4_num}")
    print(f"            = {numerator}")

    # Printing the denominator calculation
    print(f"\nDenominator = 2^({n+1}) * (1 + 3^{n})")
    print(f"            = {term1_den} * ({term2_den})")
    print(f"            = {denominator}")

    # Final result
    result_float = numerator / denominator
    print("\n--- Final Result ---")
    print(f"||T||_1 = {numerator} / {denominator}")
    print(f"||T||_1 ≈ {result_float}")
    print("="*50)

if __name__ == '__main__':
    calculate_norm_for_Jn()