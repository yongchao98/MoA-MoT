def get_prime_factorization_exponents(n):
    """
    Calculates the prime factorization of a number and returns a dictionary
    of {prime: exponent}.
    For example, for n=12, it returns {2: 2, 3: 1}.
    """
    factors = {}
    divisor = 2
    temp_n = n
    while divisor * divisor <= temp_n:
        while temp_n % divisor == 0:
            factors[divisor] = factors.get(divisor, 0) + 1
            temp_n //= divisor
        divisor += 1
    if temp_n > 1:
        factors[temp_n] = factors.get(temp_n, 0) + 1
    return factors

def solve_max_independent_events():
    """
    Calculates the largest possible value of m based on the problem's parameters.
    """
    num_dice = 100
    num_sides = 6

    # The size of the sample space is num_sides ^ num_dice.
    # The maximum number of independent events is the sum of the exponents
    # in the prime factorization of the sample space size.
    # This is equivalent to num_dice * (sum of exponents in prime factorization of num_sides).
    
    # 1. Get prime factorization of the number of sides.
    side_factors = get_prime_factorization(num_sides)
    
    # 2. Calculate the exponents for the total sample space size N = 6^100.
    # The exponent for each prime p in N is (exponent of p in 6) * 100.
    total_exponents = []
    for prime in sorted(side_factors.keys()):
        exponent_for_N = side_factors[prime] * num_dice
        total_exponents.append(exponent_for_N)
        
    # 3. The maximum number of events 'm' is the sum of these exponents.
    m = sum(total_exponents)

    print(f"The number of sides on a die is {num_sides}, and the number of dice is {num_dice}.")
    print(f"The size of the sample space is {num_sides}^{num_dice}.")
    print(f"The prime factorization of the number of sides ({num_sides}) is {side_factors}.")
    print(f"Therefore, the prime factorization of the sample space size is 2^{total_exponents[0]} * 3^{total_exponents[1]}.")
    print("The maximum number of mutually independent events is the sum of the exponents of this prime factorization.")
    
    equation_parts = [str(exp) for exp in total_exponents]
    equation_str = " + ".join(equation_parts)
    
    print(f"The largest possible value of m is: {equation_str} = {m}")

# Execute the solution
solve_max_independent_events()