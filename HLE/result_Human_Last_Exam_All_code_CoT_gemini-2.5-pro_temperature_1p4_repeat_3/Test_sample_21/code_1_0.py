def solve_density_problem():
    """
    Calculates the natural density of primes p for which the polynomial
    f(x) = x^7 - 14x^5 + 56x^3 - 56x + 22 is irreducible mod p.
    """
    
    # 1. Identify the Galois Group
    # The polynomial f(x) = x^7 - 14x^5 + 56x^3 - 56x + 22 is known in
    # algebraic number theory to have the cyclic group of order 7 as its
    # Galois group over the rational numbers.
    galois_group_name = "Cyclic group of order 7 (C_7)"
    group_order = 7
    polynomial_degree = 7

    print("Step 1: Identify the Galois group of the polynomial.")
    print(f"The polynomial is f(x) = x^7 - 14x^5 + 56x^3 - 56x + 22.")
    print(f"Its Galois group G over the rational numbers is the {galois_group_name}.")
    print(f"The order of the group, |G|, is {group_order}.")
    print("-" * 20)

    # 2. State the condition for irreducibility mod p
    print("Step 2: Determine the condition for irreducibility modulo p.")
    print(f"For a polynomial of prime degree {polynomial_degree}, it remains irreducible modulo a prime p")
    print(f"if and only if the Frobenius element Frob_p acts as a {polynomial_degree}-cycle on the roots.")
    print("-" * 20)
    
    # 3. Count the number of 7-cycles in the Galois group G = C_7
    # In a cyclic group of prime order q, any non-identity element is a generator
    # and has order q. In a transitive permutation group of degree q, an element
    # of order q must be a q-cycle.
    num_7_cycles = group_order - 1

    print("Step 3: Count the number of elements in G that are 7-cycles.")
    print(f"The group G = C_7 has {group_order} elements. The elements of order 7 are the 7-cycles.")
    print("In C_7, all elements except the identity have order 7.")
    print(f"So, the number of 7-cycles in G is {group_order} - 1 = {num_7_cycles}.")
    print("-" * 20)

    # 4. Apply the Chebotarev Density Theorem
    numerator = num_7_cycles
    denominator = group_order

    print("Step 4: Apply the Chebotarev Density Theorem.")
    print("The theorem states that the density is the ratio of the number of 7-cycles")
    print("to the total number of elements in the group G.")
    print("\nThe final equation is:")
    print(f"Density = (Number of 7-cycles in G) / |G|")
    print(f"Density = {numerator} / {denominator}")


solve_density_problem()
<<<6/7>>>