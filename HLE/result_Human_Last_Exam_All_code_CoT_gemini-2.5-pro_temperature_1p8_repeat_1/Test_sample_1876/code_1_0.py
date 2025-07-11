import math

def solve_smallest_n():
    """
    This function determines the smallest N by constructing a set of 16 distinct integers
    that can form a 4x4 multiplicative magic square.
    """
    # Exponent sets for primes 2, 3, and 5
    exponents_p2 = [0, 1, 2, 3]
    exponents_p3 = [0, 1]
    exponents_p5 = [0, 1]

    # Generate the 16 distinct numbers
    numbers = []
    for e2 in exponents_p2:
        for e3 in exponents_p3:
            for e5 in exponents_p5:
                num = (2**e2) * (3**e3) * (5**e5)
                numbers.append(num)

    # The value of N is the maximum number in this set.
    N = max(numbers)
    
    # To satisfy the prompt "output each number in the final equation",
    # we show the calculation for N.
    max_e2 = max(exponents_p2)
    max_e3 = max(exponents_p3)
    max_e5 = max(exponents_p5)
    
    print(f"The set of 16 distinct numbers is:\n{sorted(numbers)}")
    print("\nThe maximum number in this set determines the smallest N.")
    print(f"The calculation for N is based on the maximum exponents used:")
    print(f"N = 2^{max_e2} * 3^{max_e3} * 5^{max_e5} = {2**max_e2} * {3**max_e3} * {5**max_e5} = {N}")

solve_smallest_n()