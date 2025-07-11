def solve_cardinality_problem():
    """
    This function determines which of the given sets have the same cardinality as [0, 1],
    which is the cardinality of the continuum (c).

    The analysis for each set is provided in the comments. Cardinalities are represented
    by strings for clarity:
    - "aleph_0": Countably infinite (like N, Q, Z).
    - "continuum": Cardinality of the continuum, c = 2^aleph_0 (like R, [0, 1]).
    - "beth_2": Cardinality of the power set of R, 2^c.
    """
    target_cardinality = "continuum"
    correct_options = []
    
    # A. (0, 1): The open interval has a bijection to R, e.g., f(x)=tan(pi*(x-0.5)).
    # Its cardinality is the continuum.
    if "continuum" == target_cardinality:
        correct_options.append("A")

    # B. N: The set of natural numbers is the definition of countably infinite.
    # Its cardinality is aleph_0.
    if "aleph_0" == target_cardinality:
        correct_options.append("B")

    # C. Q: The set of rational numbers is countably infinite.
    # Its cardinality is aleph_0.
    if "aleph_0" == target_cardinality:
        correct_options.append("C")

    # D. R: The set of real numbers is the canonical example of a set with the
    # cardinality of the continuum.
    if "continuum" == target_cardinality:
        correct_options.append("D")

    # E. R \ Q: The set of irrational numbers. Since R = Q U (R \ Q),
    # we have c = aleph_0 + |R \ Q|. This implies |R \ Q| = c.
    if "continuum" == target_cardinality:
        correct_options.append("E")

    # F. C (Complex numbers): C is equivalent to R^2. |R^2| = c * c = c.
    if "continuum" == target_cardinality:
        correct_options.append("F")

    # G. H (Quaternions): H is equivalent to R^4. |R^4| = c^4 = c.
    if "continuum" == target_cardinality:
        correct_options.append("G")

    # H. {x: c'(x) = 0}, where c(x) is the Cantor function. This set is [0, 1]
    # without the Cantor set C. Since |C| = c and |[0, 1]| = c, the
    # cardinality of the remaining set is also c.
    if "continuum" == target_cardinality:
        correct_options.append("H")

    # I. The set of finite strings over any finite/countable alphabet is a
    # countable union of finite sets, which is countable (aleph_0).
    if "aleph_0" == target_cardinality:
        correct_options.append("I")

    # J. Set of all points in a countably infinite dimensional space (R^N).
    # Its cardinality is |R|^|N| = c^aleph_0 = (2^aleph_0)^aleph_0 = 2^(aleph_0*aleph_0) = 2^aleph_0 = c.
    if "continuum" == target_cardinality:
        correct_options.append("J")

    # K. Set of all lattice points in a countably infinite dimensional space (Z^N).
    # Its cardinality is |Z|^|N| = aleph_0^aleph_0. Since 2^aleph_0 <= aleph_0^aleph_0 <= c^aleph_0 = c, this is c.
    if "continuum" == target_cardinality:
        correct_options.append("K")

    # L. N x N: The Cartesian product of two countable sets is countable (aleph_0).
    if "aleph_0" == target_cardinality:
        correct_options.append("L")

    # M. R x R: Same as C (complex numbers). Cardinality is c * c = c.
    if "continuum" == target_cardinality:
        correct_options.append("M")

    # N. 2^N: The power set of N. By definition, its cardinality is 2^|N| = 2^aleph_0 = c.
    if "continuum" == target_cardinality:
        correct_options.append("N")

    # O. 2^Q: The power set of Q. Since |Q| = aleph_0, |2^Q| = 2^aleph_0 = c.
    if "continuum" == target_cardinality:
        correct_options.append("O")

    # P. 2^C: The power set of C. |C|=c, so |2^C| = 2^c. By Cantor's theorem, 2^c > c.
    if "beth_2" == target_cardinality:
        correct_options.append("P")

    # Q. 2^(2^Z): The power set of the power set of Z. |Z|=aleph_0, |2^Z|=c, so |2^(2^Z)|=2^c.
    if "beth_2" == target_cardinality:
        correct_options.append("Q")

    correct_options.sort()
    final_answer = "".join(correct_options)

    if not final_answer:
        print("0")
    else:
        print(final_answer)

solve_cardinality_problem()