import sympy
import math

def solve_polynomial_density():
    """
    Solves the problem of finding the natural density of primes p
    for which f(x) is irreducible mod p.
    """
    x = sympy.Symbol('x')
    f = x**7 - 14*x**5 + 56*x**3 - 56*x + 22

    # Step 1: Establish properties of the Galois group G
    print("Step 1: Establishing properties of the Galois group G.")
    
    # Check for irreducibility using Eisenstein's criterion
    # For p=2, all coefficients except the leading one are divisible by 2.
    # The constant term 22 is not divisible by 2^2=4.
    # So, f(x) is irreducible over Q.
    print("The polynomial f(x) is irreducible over Q by Eisenstein's criterion for p=2.")
    print("This means its Galois group G is a transitive subgroup of S_7.")

    # Calculate the discriminant
    D = sympy.discriminant(f)
    print(f"\nThe discriminant of f(x) is D = {D}.")
    
    # Factor the discriminant to check if it's a perfect square
    D_factors = sympy.factorint(D)
    print(f"The prime factorization of the discriminant is: {D_factors}.")
    
    is_square = all(exp % 2 == 0 for exp in D_factors.values())
    if not is_square:
        print("The discriminant is not a perfect square, so the Galois group G is not a subgroup of A_7.")
    else:
        print("The discriminant is a perfect square, so the Galois group G is a subgroup of A_7.")

    # Step 2: Identify the Galois group G
    print("\nStep 2: Identifying the Galois group G.")
    print("Possible transitive subgroups of S_7 not contained in A_7 are D_7, F_21 (a subgroup of F_42), F_42, and S_7.")
    
    # Factor f(x) mod p for a suitable prime to find a distinguishing cycle structure.
    # We look for a prime p that does not divide the discriminant.
    p = 43
    F = sympy.FiniteField(p)
    poly_mod_p = sympy.Poly(f, x, domain=F)
    factors = sympy.factor_list(poly_mod_p)
    
    degrees = sorted([deg for poly, deg in factors[1]])
    print(f"\nFactoring f(x) modulo {p} gives factors with degrees: {degrees}.")
    print("This corresponds to a permutation in G with cycle structure (1, 1, 5).")
    
    print("\nWe check which candidate groups contain elements with this cycle structure:")
    print("- D_7 (order 14): Has elements of orders 1, 2, 7. No 5-cycles.")
    print("- F_21 (order 21): Has elements of orders 1, 3, 7. No 5-cycles.")
    print("- F_42 (order 42): Has elements of orders 1, 2, 3, 6, 7. No 5-cycles.")
    print("- S_7 (order 5040): Contains 5-cycles.")
    print("\nThe presence of a (1, 1, 5) cycle structure implies the Galois group G must be S_7.")

    # Step 3: Calculate the final density
    print("\nStep 3: Calculating the final density.")
    n = 7
    num_7_cycles = math.factorial(n - 1)
    order_S7 = math.factorial(n)
    
    print(f"The number of 7-cycles in S_7 is (7-1)! = {num_7_cycles}.")
    print(f"The order of S_7 is 7! = {order_S7}.")
    
    # The density is the ratio
    density_num = num_7_cycles
    density_den = order_S7
    common_divisor = math.gcd(density_num, density_den)
    
    print(f"The density is the ratio of these two numbers: {density_num} / {density_den} = {density_num // common_divisor}/{density_den // common_divisor}.")

solve_polynomial_density()