import sympy

def analyze_curves():
    """
    Analyzes each curve to determine if it has good reduction for all primes p > 2.
    """
    x = sympy.Symbol('x')
    curves = {
        'A': sympy.Poly(x**5 + 3, x, domain='ZZ'),
        'B': sympy.Poly(x**5 - 1, x, domain='ZZ'),
        'C': sympy.Poly(x**6 - 1, x, domain='ZZ'),
        'D': sympy.Poly(2*x**5 + 2*x**3 + 1, x, domain='ZZ'),
        'E': sympy.Poly(4*x**5 + 4*x**3 + x**2 + 4*x, x, domain='ZZ')
    }

    print("--- Analysis of Each Curve ---")

    for label, f in curves.items():
        print(f"\nCurve {label}: z^2 = {sympy.pretty(f.as_expr())}")

        # Step 1: Calculate discriminant and find candidate primes
        disc = sympy.discriminant(f)
        factors = sympy.factorint(disc)
        
        print(f"  Discriminant: {disc} = {sympy.pretty(sympy.factorint(disc, visual=True))}")
        
        odd_prime_factors = [p for p in factors if p not in [-1, 2]]

        if not odd_prime_factors:
            print("  This curve has no odd prime factors in its discriminant.")
            print("  Conclusion: Curve has good reduction for all primes p > 2.")
            continue # Move to next curve

        print(f"  Candidate primes for bad reduction (p > 2): {odd_prime_factors}")

        # Step 2: Verify bad reduction using the GCD test
        has_bad_reduction_above_2 = False
        for p in odd_prime_factors:
            try:
                # Reduce polynomial and its derivative modulo p
                f_mod_p = f.set_modulus(p)
                df_mod_p = sympy.diff(f_mod_p, x)

                # Check if degree drops, which can complicate the simple discriminant test
                if f.degree() != f_mod_p.degree() or sympy.degree(f.diff(x)) != df_mod_p.degree():
                    print(f"  Note: Degree of polynomial or its derivative drops at p={p}.")

                # The definitive test: gcd of the reduced polynomial and its derivative
                the_gcd = sympy.gcd(f_mod_p, df_mod_p)
                
                print(f"  Testing at p={p}: gcd(f mod p, f' mod p) = {sympy.pretty(the_gcd.as_expr())}")

                if not the_gcd.is_one:
                    print(f"  Conclusion for p={p}: The gcd is not 1, so the curve has BAD reduction at p={p}.")
                    has_bad_reduction_above_2 = True
                    break # Found a prime of bad reduction, no need to check others for this curve
                else:
                    # This is the interesting case for curve D
                    print(f"  Conclusion for p={p}: The gcd is 1, so the curve has GOOD reduction at p={p}.")

            except Exception as e:
                print(f"  An error occurred during analysis for p={p}: {e}")

        if has_bad_reduction_above_2:
            print("  Final Conclusion: This curve does NOT have good reduction for all primes p > 2.")
        else:
            print("  Final Conclusion: This curve HAS good reduction for all primes p > 2.")


if __name__ == '__main__':
    analyze_curves()