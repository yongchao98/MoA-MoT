import math

def solve_cubic_equation():
    """
    Numerically solves the cubic equation p^3 + p - 1 = 0 for the real root in [0, 1]
    using the bisection method.
    """
    f = lambda p: p**3 + p - 1
    low, high = 0.0, 1.0
    # Iterate for a high number of times to ensure precision
    for _ in range(100):
        mid = (low + high) / 2
        if f(mid) < 0:
            low = mid
        else:
            high = mid
    return (low + high) / 2

def main():
    """
    Calculates the minimal cost coefficient and prints the process.
    """
    print("To find the minimal cost, we first need to find the optimal probability 'p' for a 'yes' answer.")
    print("This probability is the solution to the equation derived from minimizing the cost per bit of information.")

    p_cubed_coeff = 1
    p_coeff = 1
    constant_term = -1
    print(f"\nThe equation is: ({p_cubed_coeff})*p^3 + ({p_coeff})*p + ({constant_term}) = 0")
    
    # 1. Solve for the optimal probability p (alpha)
    optimal_p = solve_cubic_equation()
    print(f"\nThe real root of this equation is p = {optimal_p:.6f}")
    
    # 2. Calculate the minimal cost using the formula
    print("\nThe minimal cost per bit 'C' is given by the formula: C = -1 / log2(p)")
    
    numerator = -1
    log_val = math.log2(optimal_p)
    
    print(f"\nBreaking down the calculation:")
    print(f"C = {numerator} / log2({optimal_p:.6f})")
    print(f"C = {numerator} / {log_val:.6f}")

    cost = numerator / log_val
    print(f"\nThus, the minimal cost coefficient is {cost:.3f}")

if __name__ == "__main__":
    main()