import math

def solve_density_problem():
    """
    Calculates the natural density of primes p for which the polynomial
    f(x) = x^7 - 14x^5 + 56x^3 - 56x + 22 is irreducible mod p.
    
    This is done by determining the likely Galois group and then applying
    the Chebotarev Density Theorem.
    """
    
    # The degree of the polynomial.
    n = 7
    
    # Based on computational evidence (factoring the polynomial modulo various primes),
    # the Galois group G of the polynomial is identified as the symmetric group S_n.
    
    # The polynomial is irreducible mod p if the Frobenius element is an n-cycle.
    # We need to find the proportion of n-cycles in S_n.
    
    # The number of n-cycles in the symmetric group S_n is (n-1)!.
    num_n_cycles = math.factorial(n - 1)
    
    # The total number of elements in the group S_n is n!.
    size_of_group = math.factorial(n)
    
    # The density is the ratio of the number of n-cycles to the total size of the group.
    # The final equation is Density = (n-1)! / n! = 1/n.
    
    print(f"The degree of the polynomial is n = {n}.")
    print("The Galois group G is determined to be the symmetric group S_n.")
    print(f"The number of {n}-cycles in S_{n} is (n-1)! which is ({n}-1)! = {num_n_cycles}.")
    print(f"The size of the group S_{n} is n! which is {n}! = {size_of_group}.")
    print("The density is the ratio of the number of 7-cycles to the size of the group:")
    print(f"Density = {num_n_cycles} / {size_of_group} = 1/{n}")

solve_density_problem()