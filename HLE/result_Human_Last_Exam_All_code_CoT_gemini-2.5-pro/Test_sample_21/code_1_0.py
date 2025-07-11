import math

def calculate_density():
    """
    Calculates the natural density of primes p for which a polynomial with a
    cyclic Galois group of prime order n remains irreducible mod p.
    """
    # The degree of the polynomial and the order of the cyclic Galois group C_n.
    n = 7

    # The order of the Galois group G = C_n.
    group_order = n

    # The number of elements in G that correspond to an irreducible polynomial
    # mod p are the n-cycles. In C_n, these are the generators.
    # The number of generators of C_n is given by Euler's totient function phi(n).
    # Since n=7 is prime, phi(7) = 7 - 1 = 6.
    num_n_cycles = n - 1

    # The density is the ratio of the number of n-cycles to the group order.
    density_numerator = num_n_cycles
    density_denominator = group_order

    print("The problem is to find the density of primes p for which the polynomial is irreducible mod p.")
    print("This density is given by the ratio of the number of n-cycles in the Galois group to the order of the group.")
    print(f"For the given polynomial, the Galois group is C_{n}.")
    print(f"The order of the group is {group_order}.")
    print(f"The number of {n}-cycles (generators) in the group is phi({n}) = {num_n_cycles}.")
    print("The final equation for the density is:")
    print(f"{density_numerator} / {density_denominator}")

calculate_density()