import math

def calculate_M_z_1(n):
    """
    Calculates the approximate value of M_z(1, n) based on the derived formula.
    M_z(1, n) = (-1)^n * (2*n * comb(2*n, n) / n^n) * (8 / pi^2)^n
    """
    if n == 0:
        return float('inf')
    
    # Python's math.comb(n, k) calculates the binomial coefficient "n choose k"
    try:
        binom_coeff = math.comb(2 * n, n)
    except ValueError:
        return float('inf') # Handles cases where comb might fail for non-integer or negative n

    term1 = (-1)**n
    term2 = (2 * n * binom_coeff) / (n**n)
    term3 = (8 / (math.pi**2))**n
    
    return term1 * term2 * term3

def find_min_magnetization():
    """
    Finds the n that minimizes M_z(1, n) and calculates the value.
    """
    min_M = float('inf')
    n_min = 0

    # We test integer values of n as 'number of spins' implies n is an integer.
    # A search up to n=20 is sufficient to observe the trend.
    for n in range(1, 21):
        M_val = calculate_M_z_1(n)
        # print(f"n = {n}, M_z(1) = {M_val:.4f}") # Uncomment for debugging
        if M_val < min_M:
            min_M = M_val
            n_min = n
            
    print(f"The minimum magnetization occurs at n_min = {n_min}.")

    # For n_min = 1, the formula simplifies significantly:
    # M_z(1, 1) = (-1)^1 * (2*1 * comb(2,1) / 1^1) * (8/pi^2)^1
    #           = -1 * (2 * 2 / 1) * (8/pi^2)
    #           = -4 * 8 / pi^2
    #           = -32 / pi^2
    numerator = -32
    denominator_str = "pi^2"
    denominator_val = math.pi**2
    
    final_value = numerator / denominator_val

    print(f"The value is given by the equation: {numerator} / {denominator_str}")
    print(f"Numerical value: {final_value}")

find_min_magnetization()

<<< -3.242364041131102