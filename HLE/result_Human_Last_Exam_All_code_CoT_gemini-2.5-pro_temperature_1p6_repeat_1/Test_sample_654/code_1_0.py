import math

def solve_polynomial_coefficients():
    """
    This function calculates the number of coefficients not divisible by p^k
    in the given recursively defined polynomial problem.

    The problem asks for the number of coefficients not divisible by p^k in the polynomial
    Q(x) = F_{p, k}(...F_{p, k}(x)...), where the function F is applied p^n times.

    The derivation shows that the final answer is given by the formula:
    Result = p^min(n, k-1) + 1
    """

    # User-defined values for p, k, n
    # According to the problem statement:
    # p is an odd prime (>= 3)
    # k and n are integers (>= 1)
    
    # You can change these values to test other cases.
    p = 3
    k = 1
    n = 1
    
    print(f"Given parameters are p = {p}, k = {k}, n = {n}.")
    print("The formula for the number of coefficients not divisible by p^k is:")
    print("Result = p^(min(n, k-1)) + 1")
    
    # Step 1: Calculate min(n, k-1)
    min_val = min(n, k - 1)
    print(f"\nFirst, we calculate the exponent, which is min(n, k-1):")
    print(f"min({n}, {k}-1) = min({n}, {k-1}) = {min_val}")

    # Step 2: Calculate p to the power of min_val
    power_val = int(math.pow(p, min_val))
    print(f"\nNext, we calculate p to the power of this exponent:")
    print(f"{p}^{min_val} = {power_val}")

    # Step 3: Add 1 to get the final result
    result = power_val + 1
    print(f"\nFinally, we add 1 to get the result:")
    print(f"{power_val} + 1 = {result}")

    print("\n---")
    print("The final equation with numbers is:")
    print(f"{p}^(min({n}, {k-1})) + 1 = {p}^{min_val} + 1 = {power_val} + 1 = {result}")

if __name__ == '__main__':
    solve_polynomial_coefficients()
