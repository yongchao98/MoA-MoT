import math

def get_prime_factorization(n):
    """
    Returns the prime factorization of n as a dictionary of {prime: exponent}.
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

def solve_dice_problem():
    """
    Calculates the largest possible number of mutually independent events for
    rolling a number of dice.
    """
    num_dice = 100
    sides_per_die = 6

    print("This script calculates the largest possible value of m (number of mutually independent events).")
    print("The experiment is rolling {} regular {}-sided dice.".format(num_dice, sides_per_die))
    print("-" * 30)

    # Explain the theory
    print("The size of the sample space is given by (sides_per_die)^num_dice.")
    print("In our case, this is {}^({}).".format(sides_per_die, num_dice))
    print("\nA theorem in probability states that the maximum number of mutually independent events")
    print("on a finite sample space is the sum of the exponents in the prime factorization of its size.")
    print("-" * 30)

    # Step 1: Factorize the number of sides
    print("Step 1: Find the prime factorization of the number of sides on a single die.")
    side_factors = get_prime_factorization(sides_per_die)
    side_factor_str = " * ".join([f"{p}^{e}".replace("^1", "") for p, e in side_factors.items()])
    print(f"The prime factorization of {sides_per_die} is: {side_factor_str}")
    print("-" * 30)

    # Step 2: Determine exponents for the total sample space
    print("Step 2: Determine the exponents for the prime factorization of the total sample space size, {}^({}).".format(sides_per_die, num_dice))
    print(f"Since Size = ({side_factor_str})^{num_dice}, the exponents of the primes are multiplied by {num_dice}.")
    
    total_space_exponents = {p: e * num_dice for p, e in side_factors.items()}
    total_space_exponents_str = " * ".join([f"{p}^{e}" for p, e in total_space_exponents.items()])
    print(f"The prime factorization of {sides_per_die}^{num_dice} is: {total_space_exponents_str}")
    print("-" * 30)

    # Step 3: Sum the exponents
    print("Step 3: Sum the exponents to find the maximum number of events (m).")
    exponents = list(total_space_exponents.values())
    max_m = sum(exponents)

    equation_parts = [str(e) for e in exponents]
    print("The final calculation is:")
    print(f"{' + '.join(equation_parts)} = {max_m}")
    print("-" * 30)
    print(f"The largest possible value of m is {max_m}.")


if __name__ == "__main__":
    solve_dice_problem()