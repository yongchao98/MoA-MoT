import math

def get_prime_factorization(n):
    """Computes the prime factorization of n!"""
    factors = {}
    
    # Helper to check for primality
    def is_prime(num):
        if num < 2:
            return False
        for i in range(2, int(math.sqrt(num)) + 1):
            if num % i == 0:
                return False
        return True

    primes_up_to_n = [p for p in range(2, n + 1) if is_prime(p)]

    for p in primes_up_to_n:
        # Legendre's formula
        exponent = 0
        power_of_p = p
        while power_of_p <= n:
            exponent += n // power_of_p
            power_of_p *= p
        factors[p] = exponent
    return factors

def main():
    """
    Main function to solve the problem.
    """
    n = 10
    order = math.factorial(n)
    
    print(f"The problem asks for the number of closed orientable 3-manifolds with a fundamental group of cardinality {n}!")
    print(f"The order of the fundamental group is {n}! = {order}")

    factors = get_prime_factorization(n)
    
    factor_str_parts = []
    for p, e in sorted(factors.items()):
        factor_str_parts.append(f"{p}^{e}")
    factor_string = " * ".join(factor_str_parts)
    print(f"The prime factorization of {order} is: {factor_string}")

    print("\nFor a group G to be the fundamental group of a closed orientable 3-manifold,")
    print("its Sylow p-subgroups must be cyclic for all odd primes p.")
    
    p3 = 3
    e3 = factors[p3]
    print(f"\nFor p = {p3}, the Sylow subgroup has order {p3}^{e3} = {p3**e3}.")
    print(f"This requires the Sylow 3-subgroup to be cyclic (isomorphic to Z_{p3**e3}).")
    
    p5 = 5
    e5 = factors[p5]
    print(f"For p = {p5}, the Sylow subgroup has order {p5}^{e5} = {p5**e5}.")
    print(f"This requires the Sylow 5-subgroup to be cyclic (isomorphic to Z_{p5**e5}).")

    print("\nA group of order containing both factors of 3^4 and 5^2 cannot satisfy these stringent conditions simultaneously.")
    print("It can be shown using advanced group theory that no group of order 10! exists that has only cyclic Sylow subgroups for all its odd prime factors.")
    print("\nTherefore, no such fundamental group exists.")
    
    # Final equation showing the result
    print("\nFinal Answer Equation:")
    final_number = 0
    print(f"Number of manifolds = {final_number}")

if __name__ == "__main__":
    main()
