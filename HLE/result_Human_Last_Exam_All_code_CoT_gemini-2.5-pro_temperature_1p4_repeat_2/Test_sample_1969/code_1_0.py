import sys

def get_fair_division_formula(n):
    """
    Calculates and prints the formula for the Shapley value c_k for a given n.

    Args:
        n (int): The number of people in the band.
    """
    if not isinstance(n, int) or n <= 1:
        print("Please provide an integer n > 1.", file=sys.stderr)
        return

    # Calculate terms that depend only on n.
    # Note: These calculations result in integers because n(n+1) is always even,
    # and n^2(n+1)^2 is divisible by 4.
    term1 = n**2 * (n+1)**2 // 24
    term2 = (3*n + 1) * (n+2)
    
    # Coefficients for the quadratic formula c_k = A*k - B*k^2
    # c_k = (term1 * term2) * k - (term1 * 6) * k^2
    coeff_k = term1 * term2
    coeff_k2 = term1 * 6

    # Output the formula for c_k.
    # The problem asks to "output each number in the final equation".
    print(f"For n = {n}, the fair amount for person p_k is given by the formula:")
    print(f"c_k = {coeff_k}*k - {coeff_k2}*k^2")
    print("\nFor example:")
    
    # Calculate for k=1 and k=n to show some concrete values
    c_1 = coeff_k * 1 - coeff_k2 * 1**2
    c_n = coeff_k * n - coeff_k2 * n**2
    print(f"Person p_1 receives: c_1 = {c_1:,.2f} dollars")
    print(f"Person p_{n} receives: c_{n} = {c_n:,.2f} dollars")


if __name__ == '__main__':
    # You can change the value of n here to test.
    # The problem has a constraint n > 1.
    try:
        # Example for n=10. For the interactive prompt, you can imagine the user provides this value.
        n_value = 10
        get_fair_division_formula(n_value)
    except (ValueError, IndexError):
        print("Usage: please provide an integer n > 1.", file=sys.stderr)
