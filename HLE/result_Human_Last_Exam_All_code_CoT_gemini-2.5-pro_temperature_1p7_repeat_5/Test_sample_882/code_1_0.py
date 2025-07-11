import math

def solve_for_p():
    """
    Solves the equilibrium equation for p and calculates the final result.
    """
    N = 8
    M = 3 * N

    print("The equilibrium equation for p is derived from Pi_D(p) = Pi_C(p).")
    print(f"For N={N}, M={M}, the equation is approximately:")
    print(f"(1 - p/{N})^{M} + p - 1 = 0")
    print(f"(1 - p/8)^24 + p - 1 = 0\n")

    # Define the function for which we want to find the root
    def f(p):
        return (1 - p / N)**M + p - 1

    # Bisection method to find the root
    # We know a root exists between 0 and 1.
    # f(0) = 0, but the derivative f'(0) is negative, so we search away from 0.
    # By testing, we find f(0.5) < 0 and f(1) > 0.
    low = 0.5
    high = 1.0
    
    # Iterate to find p with high precision
    for _ in range(100):
        mid = (low + high) / 2
        if f(mid) < 0:
            low = mid
        else:
            high = mid
            
    p = (low + high) / 2
    
    print(f"The calculated probability p is approximately: {p:.6f}")
    
    one_minus_p = 1 - p
    print(f"The value of (1-p) is approximately: {one_minus_p:.6f}")
    
    result = 10000 * one_minus_p
    print(f"The value of 10000 * (1-p) is: {result:.6f}")
    
    final_answer = math.floor(result)
    print(f"The floor of this value is: {final_answer}")
    
    return final_answer

final_value = solve_for_p()
# The problem asks for the answer in a specific format at the end.
# The integer value is what needs to be formatted.
print(f'<<<{final_value}>>>')
