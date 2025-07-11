def solve_chromatic_roots_problem():
    """
    Analyzes statements about chromatic and orbital chromatic roots and identifies the true ones.

    The analysis for each statement is as follows:

    A. Real orbital chromatic roots are bounded by the greatest real chromatic root.
       - This is a known theorem by Cameron and Jackson (2007). They proved that the real
         roots of the orbital chromatic polynomial P_{G,Γ}(x) are indeed bounded above
         by the largest real root of the standard chromatic polynomial P_G(x).
       - Status: True.

    B. Chromatic roots may not be real.
       - Chromatic polynomials have integer coefficients. The fundamental theorem of algebra
         guarantees that their roots are in the complex plane, but not necessarily on the
         real line. In fact, many graphs are known to have non-real complex chromatic roots.
         A strong result by Sokal (2004) showed that chromatic roots are dense in the
         entire complex plane.
       - Status: True.

    C. Real chromatic roots may take on negative values.
       - It is a well-known property that the coefficients of the chromatic polynomial P_G(x),
         when written in the standard basis {x^k}, alternate in sign.
         P_G(x) = x^n - a_{n-1}x^{n-1} + a_{n-2}x^{n-2} - ... where a_i >= 0.
         This structure implies that for any x < 0, all terms in the polynomial are positive
         (or all negative, depending on n), so P_G(x) cannot be zero. Therefore, there are
         no negative real chromatic roots.
       - Status: False.

    D. Real chromatic roots may take on non-integer values.
       - While many simple graphs have integer roots, this is not true in general. A famous
         result by W. T. Tutte showed that for a specific planar triangulation, the value
         (3 + sqrt(5))/2 ≈ 2.618 is a chromatic root. This value is real but not an integer.
       - Status: True.

    E. Chromatic polynomials may have roots between 0 and 1 (i.e., r ∈ (0, 1)).
       - This is known to be false. A significant result in the theory, first proved by
         Jackson (1993) and later given a simpler proof by Thomassen (1997), states that
         chromatic polynomials have no roots in the interval (0, 1).
       - Status: False.
    """
    
    # Store the boolean result for each statement
    results = {
        'A': True,
        'B': True,
        'C': False,
        'D': True,
        'E': False
    }

    # Filter for true statements
    true_statements = [key for key, value in results.items() if value]

    # Sort the letters alphabetically
    true_statements.sort()

    # Join into a single string
    final_answer = "".join(true_statements)

    print(final_answer)

solve_chromatic_roots_problem()