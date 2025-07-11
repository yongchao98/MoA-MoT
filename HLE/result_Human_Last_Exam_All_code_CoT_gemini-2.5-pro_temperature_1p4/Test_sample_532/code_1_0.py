import math

def get_prime_factorization(n):
    """
    Calculates the prime factorization of a positive integer n.
    Returns a dictionary where keys are prime factors and values are their powers.
    """
    factors = {}
    d = 2
    temp = n
    while d * d <= temp:
        while (temp % d) == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

def classify_nilpotent_group_by_order(order):
    """
    Classifies a finite nilpotent group as filled or not based on its order.
    A finite nilpotent group is filled if and only if it is a p-group,
    which means its order is a power of a single prime.
    """
    print(f"--- Analyzing a nilpotent group of order {order} ---")
    if not isinstance(order, int) or order < 1:
        print("Error: The order must be a positive integer.")
        return

    if order == 1:
        print("The group is the trivial group {e}.")
        print("The order is 1, which can be written as p^0 for any prime p.")
        print("Final Equation: 1 = 2^0")
        print("Result: The trivial group is a p-group and therefore is a filled nilpotent group.")
        return

    factors = get_prime_factorization(order)
    
    # "Remember in the final code you still need to output each number in the final equation!"
    # This is interpreted as showing the prime factorization of the order.
    factor_list = [f"{p}^{e}" for p, e in sorted(factors.items())]
    final_equation = f"{order} = {' * '.join(factor_list)}"
    print(f"Final Equation: {final_equation}")

    if len(factors) == 1:
        p, a = list(factors.items())[0]
        print(f"The order is a power of a single prime ({p}).")
        print(f"The group is a p-group.")
        print("Result: The group is a filled nilpotent group.")
    else:
        primes = sorted(factors.keys())
        print(f"The order is divisible by multiple distinct primes: {', '.join(map(str, primes))}.")
        print(f"The group is not a p-group and has a quotient of the form C_{primes[0]} x C_{primes[1]}.")
        print("Result: The group is NOT a filled nilpotent group.")

# --- Example Usage ---
# A group of order 32 = 2^5 is a 2-group, hence filled.
classify_nilpotent_group_by_order(32)
print("\n")
# A group of order 12 = 2^2 * 3 is not a p-group, hence not filled.
# For example, C_12 is nilpotent but not filled.
classify_nilpotent_group_by_order(12)
print("\n")
# A group of order 243 = 3^5 is a 3-group, hence filled.
classify_nilpotent_group_by_order(243)
print("\n")
# A group of order 30 = 2 * 3 * 5 is not a p-group, hence not filled.
classify_nilpotent_group_by_order(30)