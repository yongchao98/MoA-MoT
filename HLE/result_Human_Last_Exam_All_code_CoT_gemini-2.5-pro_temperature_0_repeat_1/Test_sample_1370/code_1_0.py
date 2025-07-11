def get_prime_factorization_exponents(n):
    """
    Calculates the prime factorization of a number n and returns a dictionary
    of {prime: exponent}.
    For example, for n=12, it returns {2: 2, 3: 1}.
    """
    factors = {}
    d = 2
    temp_n = n
    while d * d <= temp_n:
        while temp_n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp_n //= d
        d += 1
    if temp_n > 1:
        factors[temp_n] = factors.get(temp_n, 0) + 1
    return factors

def solve_largest_m():
    """
    Solves the problem of finding the largest possible value of m.
    """
    num_dice = 100
    num_sides = 6

    print(f"The experiment consists of rolling {num_dice} dice, each with {num_sides} sides.")
    print(f"The total number of outcomes in the sample space is N = {num_sides}^{num_dice}.")
    print("-" * 30)

    print("To find the maximum number of mutually independent events, we use a theorem")
    print("that relates this number to the prime factorization of the sample space size N.")
    print("\nStep 1: Find the prime factorization of the number of sides of a single die.")
    
    side_factors = get_prime_factorization_exponents(num_sides)
    
    factor_parts = []
    for p, e in side_factors.items():
        factor_parts.append(f"{p}^{e}")
    print(f"The prime factorization of {num_sides} is: {' * '.join(factor_parts)}")
    
    print("\nStep 2: Determine the prime factorization of the total sample space size, N.")
    print(f"N = {num_sides}^{num_dice} = ({' * '.join(factor_parts)})^{num_dice}")

    space_factor_parts = []
    exponent_sum_parts = []
    total_exponents = 0
    for p, e in side_factors.items():
        total_exponent = e * num_dice
        space_factor_parts.append(f"{p}^{total_exponent}")
        exponent_sum_parts.append(str(total_exponent))
        total_exponents += total_exponent
        
    print(f"So, the prime factorization of N is: {' * '.join(space_factor_parts)}")

    print("\nStep 3: Apply the theorem.")
    print("The theorem states that the maximum number of mutually independent events (m)")
    print("is the sum of the exponents in the prime factorization of N.")
    
    print("\nStep 4: Calculate the final answer.")
    print("m = " + " + ".join(exponent_sum_parts))
    print(f"m = {total_exponents}")
    
    return total_exponents

if __name__ == "__main__":
    max_m = solve_largest_m()
    print(f"\n<<<The largest possible value of m is {max_m}.>>>")
