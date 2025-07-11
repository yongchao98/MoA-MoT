import collections

def solve_dice_problem():
    """
    Calculates the largest possible number of mutually independent events
    for rolling 100 six-sided dice.
    """
    base = 6
    power = 100

    # Step 1: Prime factorize the base number
    factors = collections.defaultdict(int)
    d = 2
    temp_n = base
    while d * d <= temp_n:
        while temp_n % d == 0:
            factors[d] += 1
            temp_n //= d
        d += 1
    if temp_n > 1:
        factors[temp_n] += 1

    # Step 2: Calculate the exponents for the full number N = base^power
    final_exponents = []
    for prime in sorted(factors.keys()):
        exponent = factors[prime]
        final_exponents.append(exponent * power)
        
    # Step 3: Sum the exponents and print the final equation
    total_m = sum(final_exponents)
    
    equation_parts = [str(e) for e in final_exponents]
    equation_str = " + ".join(equation_parts)
    
    print(f"{equation_str} = {total_m}")

solve_dice_problem()