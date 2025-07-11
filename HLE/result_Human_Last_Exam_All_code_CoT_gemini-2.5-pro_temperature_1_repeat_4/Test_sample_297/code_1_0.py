def solve():
    """
    Calculates the total number of multiplication operations in the fully expanded
    sum-of-products expression for the sum bit s2 of a 3-bit binary adder.
    """
    # The sum-of-products expression for c2 is: c2 = a1*b1 + a0*b0*a1'*b1 + a0*b0*a1*b1'
    # This expression has 3 product terms. We list the number of literals in each term.
    # Term 1: a1*b1 (2 literals)
    # Term 2: a0*b0*a1'*b1 (4 literals)
    # Term 3: a0*b0*a1*b1' (4 literals)
    c2_term_literals = [2, 4, 4]

    # The sum-of-products expression for c2' is:
    # c2' = a0'*a1' + a0'*b1' + b0'*a1' + b0'*b1' + a1'*b1'
    # This expression has 5 product terms, each with 2 literals.
    c2_prime_term_literals = [2, 2, 2, 2, 2]

    # The final expression for s2 is: s2 = a2'b2'c2 + a2'b2c2' + a2b2'c2' + a2b2c2
    # We calculate the number of multiplications for each of these four major components.
    # A product of n literals has (n-1) multiplications.

    # Component 1: a2' * b2' * c2
    # Expands into terms like a2'*b2'*<term_from_c2>.
    # Adds 2 literals (a2', b2') to each term of c2.
    comp1_mults = sum((2 + L - 1) for L in c2_term_literals)

    # Component 2: a2' * b2 * c2'
    # Adds 2 literals (a2', b2) to each term of c2'.
    comp2_mults = sum((2 + L - 1) for L in c2_prime_term_literals)

    # Component 3: a2 * b2' * c2'
    # Adds 2 literals (a2, b2') to each term of c2'.
    comp3_mults = sum((2 + L - 1) for L in c2_prime_term_literals)

    # Component 4: a2 * b2 * c2
    # Adds 2 literals (a2, b2) to each term of c2.
    comp4_mults = sum((2 + L - 1) for L in c2_term_literals)

    # The prompt requires printing each number in the final equation.
    # We interpret this as printing the number of multiplications contributed by each
    # of the four major components of the s2 expression.
    print(f"Number of multiplications from a2'*b2'*c2: {comp1_mults}")
    print(f"Number of multiplications from a2'*b2*c2': {comp2_mults}")
    print(f"Number of multiplications from a2*b2'*c2': {comp3_mults}")
    print(f"Number of multiplications from a2*b2*c2: {comp4_mults}")

    total_multiplications = comp1_mults + comp2_mults + comp3_mults + comp4_mults
    print(f"\nTotal number of multiplications for s2: {total_multiplications}")

solve()