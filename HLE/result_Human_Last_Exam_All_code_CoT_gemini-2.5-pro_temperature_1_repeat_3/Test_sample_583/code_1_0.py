import math

def calculate_limit_numerically():
    """
    This function calculates the value of f(n)/(n*log2(n)) for large n
    to numerically verify the analytical limit.
    """

    # Define the function f(n) based on the derived formula
    def f(n):
        if n < 1:
            return 0
        if n == 1:
            return 1
        # For n>=2, the formula is n + sum_{k=1 to n-1} floor(log2(k))
        log_sum = sum(math.floor(math.log2(k)) for k in range(1, n))
        return n + log_sum

    # Define the expression g(n) whose limit we are evaluating
    def g(n):
        if n <= 1:
            return float('nan')
        return f(n) / (n * math.log2(n))

    print("This script numerically evaluates the limit of f(n)/(n*log2(n)) as n -> infinity.")
    print("Based on analysis, the function f(n) is believed to be n + sum_{k=1 to n-1} floor(log2(k)).")
    print("-" * 30)
    print("Value of the expression for increasing n:")

    # Calculate and print g(n) for large powers of 10
    for i in range(2, 8):
        n = 10**i
        value = g(n)
        print(f"For n = {n:1.0e}, the value is: {value:.8f}")

    # State the analytical result
    limit_part1 = 0
    limit_part2 = 1
    total_limit = limit_part1 + limit_part2

    print("-" * 30)
    print("The analytical limit is calculated as follows:")
    print("lim f(n)/(n*log2(n)) = lim n/(n*log2(n)) + lim (sum_{k=1 to n-1} floor(log2(k)))/(n*log2(n))")
    print(f"The limit of the first term is {limit_part1}.")
    print(f"The limit of the second term is {limit_part2}.")
    print(f"The final limit is the sum of these two values.")
    print(f"Final Limit = {limit_part1} + {limit_part2} = {total_limit}")


calculate_limit_numerically()