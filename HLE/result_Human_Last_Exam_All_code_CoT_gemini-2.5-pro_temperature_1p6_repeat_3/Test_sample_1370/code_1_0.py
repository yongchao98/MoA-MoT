import math

def get_prime_factorization(n):
    """
    Returns a dictionary representing the prime factorization of n.
    e.g., for n=12, returns {2: 2, 3: 1}
    """
    factors = {}
    d = 2
    temp = n
    while d * d <= temp:
        while temp % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

def solve():
    """
    Solves the problem of finding the maximum number of mutually independent events.
    """
    num_dice = 100
    sides_per_die = 6

    # Step 1: Explain the theory.
    print("The problem is to find the maximum number of mutually independent events for rolling 100 dice.")
    print("This number is determined by the prime factorization of the total number of outcomes in the sample space.")
    print("-" * 20)

    # Step 2: Calculate the size of the sample space and its prime factorization.
    print(f"There are {num_dice} dice, each with {sides_per_die} sides.")
    print(f"The total number of outcomes is |Ω| = {sides_per_die}^{num_dice}.")
    
    prime_factors_of_sides = get_prime_factorization(sides_per_die)
    
    print(f"The prime factorization of the number of sides ({sides_per_die}) is: ", end="")
    factors_str = " * ".join([f"{p}^{e}" for p, e in prime_factors_of_sides.items()])
    print(factors_str)

    print(f"Therefore, the prime factorization of the sample space size |Ω| is ({factors_str})^{num_dice}.")

    total_exponents = {}
    for p, e in prime_factors_of_sides.items():
        total_exponents[p] = e * num_dice
    
    factors_str_total = " * ".join([f"{p}^{e}" for p, e in total_exponents.items()])
    print(f"|Ω| = {factors_str_total}")
    print("-" * 20)
    
    # Step 3: Sum the exponents to find the maximum number of events.
    print("The maximum number of mutually independent events (m) is the sum of the exponents in the prime factorization of |Ω|.")
    
    m = sum(total_exponents.values())
    
    exponents = list(total_exponents.values())
    
    equation_parts = [str(e) for e in exponents]
    equation_str = " + ".join(equation_parts)

    # Step 4: Output the final equation and the result.
    print(f"The exponents are {', '.join(equation_parts)}.")
    print("The final equation is:")
    print(f"{equation_str} = {m}")
    
    print("-" * 20)
    print(f"The largest possible value of m is {m}.")

solve()