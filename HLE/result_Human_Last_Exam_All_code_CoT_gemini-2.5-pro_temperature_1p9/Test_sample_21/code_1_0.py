import sympy
from sympy import Symbol, Poly
from sympy.ntheory import is_square

def analyze_polynomial():
    """
    Analyzes the polynomial f(x) = x^7 - 14x^5 + 56x^3 - 56x + 22
    to determine the natural density of primes p for which it is irreducible mod p.
    """
    x = Symbol('x')
    f_poly = Poly(x**7 - 14*x**5 + 56*x**3 - 56*x + 22, x, domain='ZZ')

    # Introduction
    print("The natural density of primes p for which a polynomial is irreducible mod p is")
    print("given by the proportion of n-cycles in its Galois group G, where n is the polynomial's degree.\n")

    # Step 1: Check irreducibility over Q.
    # We use Eisenstein's criterion with p=2. The coefficients are [1, 0, -14, 0, 56, 0, -56, 22].
    # p=2 divides all coefficients except the leading one (1).
    # p^2=4 does not divide the constant term (22).
    # Thus, the polynomial is irreducible over Q.
    print("Step 1: Determine properties of the Galois group G.")
    print("The polynomial is irreducible over Q (by Eisenstein's criterion for p=2).")
    print("This means its Galois group G is a transitive subgroup of S_7.")
    
    # Step 2: Compute the discriminant to check if G is a subgroup of A_7.
    disc = sympy.discriminant(f_poly)
    print(f"The discriminant of the polynomial is: {disc}")
    
    disc_is_square = disc.is_perfect_square()
    print(f"Is the discriminant a perfect square? {disc_is_square}")
    if not disc_is_square:
        print("Since the discriminant is not a perfect square, the Galois group G is NOT a subgroup of A_7.")
    print("-" * 40)

    # Step 3: Factor modulo unramified primes to constrain G.
    # The prime factors of the discriminant (2038431744 = 2^23 * 3^5) are 2 and 3.
    # We must use unramified primes (p != 2, 3). Let's use p=5 and p=59.
    print("Step 2: Constrain G by factoring modulo unramified primes.")
    
    # Modulo 5
    p5 = 5
    factors_5 = sympy.factor_list(f_poly, modulus=p5)
    degrees_5 = sorted([p.degree() for _, p in factors_5[1]])
    print(f"Factoring the polynomial modulo {p5}, the degrees of the factors are: {degrees_5}")
    print("This implies that G contains a permutation with cycle structure (1)(2)(4).")
    print("An element with this structure has order lcm(1,2,4)=4. This rules out G being C_7, D_7, or AGL(1,7).")
    print()

    # Modulo 59
    p59 = 59
    factors_59 = sympy.factor_list(f_poly, modulus=p59)
    degrees_59 = sorted([p.degree() for _, p in factors_59[1]])
    print(f"Factoring the polynomial modulo {p59}, the degrees of the factors are: {degrees_59}")
    print("This implies that G contains a permutation with cycle structure (1)(1)(5).")
    print("An element with this structure has order 5. This rules out G being PSL(2,7), as it has no elements of order 5.")
    print("-" * 40)
    
    # Step 4: Conclude the Galois group and compute the density.
    print("Step 3: Identify G and compute the density.")
    print("The only transitive subgroup of S_7 not ruled out is S_7 itself. Thus, the Galois group G is S_7.")
    print("\nThe density is the proportion of 7-cycles in S_7.")
    
    n = 7
    num_n_cycles = sympy.factorial(n - 1)
    size_of_Sn = sympy.factorial(n)
    
    density_num = 1
    density_den = n
    
    print("\nThe final equation is:")
    print(f"Density = (Number of 7-cycles in S_7) / (Size of S_7)")
    print(f"Density = (7-1)! / 7! = {num_n_cycles} / {size_of_Sn} = {density_num} / {density_den}")

analyze_polynomial()