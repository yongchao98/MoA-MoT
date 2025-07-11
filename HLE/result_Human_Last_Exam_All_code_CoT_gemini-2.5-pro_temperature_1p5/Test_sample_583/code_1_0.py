import math

def f(n):
    """
    Calculates the value of f(n) based on the formula from the
    solution to IMO Shortlist 2010 problem N6.
    f(n) is asymptotically equivalent to sum_{k=1 to n} floor(log2(k)).
    """
    s = 0
    # We can start from k=2 since floor(log2(1)) is 0.
    for k in range(2, n + 1):
        s += math.floor(math.log2(k))
    # The exact formula is 1 + sum_{k=1}^{n-1}, but for large n,
    # this difference is negligible for the limit calculation.
    return s

def calculate_limit_expression(n):
    """
    Calculates and prints the components and result of the expression
    f(n) / (n * log2(n)) for a given n.
    """
    # Calculate each part of the expression
    numerator_f_n = f(n)
    log_2_n = math.log2(n)
    denominator_n_log2_n = n * log_2_n
    
    # Ensure denominator is not zero
    if denominator_n_log2_n == 0:
        print(f"Cannot calculate for n = {n} as denominator is zero.")
        return

    result = numerator_f_n / denominator_n_log2_n

    # Output the numbers in the equation
    print(f"For n = {n}:")
    print(f"  f(n) = {numerator_f_n}")
    print(f"  n = {n}")
    print(f"  log2(n) = {log_2_n:.6f}")
    print(f"  The value of the expression f(n) / (n * log2(n)) is:")
    print(f"  {numerator_f_n} / ({n} * {log_2_n:.6f}) = {result:.8f}\n")


def main():
    """
    Main function to run the numerical investigation for several large values of n.
    """
    print("This program numerically investigates the limit of f(n)/(n*log2(n)) as n -> infinity.")
    print("Based on symbolic analysis, the limit should be 1.")
    print("-" * 70)
    
    # Calculate for increasingly large values of n
    for i in range(3, 8):
        n = 10**i
        calculate_limit_expression(n)
    
    print("As n increases, the value of the expression gets closer to 1.")

main()
