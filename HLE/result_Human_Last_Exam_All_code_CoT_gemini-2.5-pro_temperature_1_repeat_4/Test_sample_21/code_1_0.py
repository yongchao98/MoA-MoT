import math

def solve_polynomial_density():
    """
    This function explains and calculates the theoretical density of primes p
    for which the polynomial f(x) = x^7 - 14x^5 + 56x^3 - 56x + 22
    remains irreducible modulo p.
    """
    
    # Step 1: State the problem and the relevant theorem.
    print("The problem is to find the natural density of primes p for which")
    print("the polynomial f(x) = x^7 - 14x^5 + 56x^3 - 56x + 22 is irreducible mod p.")
    print("This density is determined by the Chebotarev Density Theorem.")
    print("-" * 20)

    # Step 2: Identify the Galois group of the polynomial.
    # The polynomial is a known example from Galois theory.
    print("The Galois group G of this polynomial over the rational numbers is the")
    print("Dihedral group D_7.")
    order_of_group = 14
    print(f"The order of D_7 is {order_of_group}.")
    print("-" * 20)
    
    # Step 3: Relate irreducibility to the group structure.
    print("For f(x) to be irreducible mod p, the corresponding Frobenius element")
    print("in the Galois group must act as a 7-cycle on the polynomial's roots.")
    print("-" * 20)
    
    # Step 4: Count the number of such elements in the group.
    # In D_7 = <r, s | r^7 = 1, s^2 = 1, srs = r^-1>, the elements that act as
    # 7-cycles are the elements of order 7. These are r, r^2, r^3, r^4, r^5, and r^6.
    num_7_cycles = 6
    print(f"The number of elements of order 7 (which are 7-cycles) in D_7 is {num_7_cycles}.")
    print("-" * 20)

    # Step 5: Calculate the final density.
    print("The density is the ratio of the number of 7-cycles to the order of the group.")
    
    # Simplify the fraction using math.gcd
    common_divisor = math.gcd(num_7_cycles, order_of_group)
    simplified_numerator = num_7_cycles // common_divisor
    simplified_denominator = order_of_group // common_divisor
    
    print("\nThe final equation for the density is:")
    print(f"{num_7_cycles} / {order_of_group} = {simplified_numerator} / {simplified_denominator}")

# Execute the function to print the explanation and result.
solve_polynomial_density()
