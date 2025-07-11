import collections

def get_prime_factorization_exponents(n):
    """
    Returns a list of exponents in the prime factorization of n.
    For example, for 12 = 2^2 * 3^1, it returns [2, 1].
    """
    factors = collections.defaultdict(int)
    d = 2
    temp = n
    while d * d <= temp:
        while (temp % d) == 0:
            factors[d] += 1
            temp //= d
        d += 1
    if temp > 1:
       factors[temp] += 1
    return list(factors.values())

def solve_dice_independence():
    """
    Calculates the maximum number of mutually independent events for rolling n dice with s sides.
    """
    num_dice = 100
    num_sides = 6

    # The size of the sample space is num_sides^num_dice.
    # The prime factorization of the size is (p1^a1 * p2^a2 * ...)^num_dice
    # which is p1^(a1*num_dice) * p2^(a2*num_dice) * ...
    # The max number of independent events is the sum of the exponents:
    # a1*num_dice + a2*num_dice + ...

    # First, get the exponents of the prime factors for a single die.
    # For 6, the prime factors are 2^1 and 3^1, so the exponents are [1, 1].
    side_exponents = get_prime_factorization_exponents(num_sides)

    # The exponents for the entire sample space are the side exponents times the number of dice.
    total_space_exponents = [exp * num_dice for exp in side_exponents]

    # The maximum number of events is the sum of these new exponents.
    max_m = sum(total_space_exponents)

    # Format the output string for the equation
    equation_parts = " + ".join(map(str, total_space_exponents))
    
    # Print each number in the final equation.
    print(f"{total_space_exponents[0]} + {total_space_exponents[1]} = {max_m}")

solve_dice_independence()
<<<200>>>