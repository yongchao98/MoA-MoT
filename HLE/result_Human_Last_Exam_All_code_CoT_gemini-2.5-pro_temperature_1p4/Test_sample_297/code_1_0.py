def solve():
    """
    Calculates the number of multiplications in the fully expanded expression for s2.
    """

    # c2 = a1*b1 + a1*a0*b0 + b1*a0*b0
    # List of the number of literals in each product term of c2's expression
    c2_term_literals = [2, 3, 3]

    # c2' = a1'*b1' + a1'*a0' + a1'*b0' + b1'*a0' + b1'*b0'
    # List of the number of literals in each product term of c2's expression
    c2_prime_term_literals = [2, 2, 2, 2, 2]

    # s2 = a2'b2'c2 + a2b2c2 + a2'b2c2' + a2b2'c2'
    # We calculate the number of multiplications for each of these four parts.
    # The prefix (e.g., a2'b2') adds 2 literals to each term from c2 or c2'.

    # Part 1: a2'b2'c2
    m1 = 0
    prefix_literals = 2
    for literals in c2_term_literals:
        total_literals = prefix_literals + literals
        m1 += (total_literals - 1)

    # Part 2: a2b2c2
    m2 = 0
    prefix_literals = 2
    for literals in c2_term_literals:
        total_literals = prefix_literals + literals
        m2 += (total_literals - 1)

    # Part 3: a2'b2c2'
    m3 = 0
    prefix_literals = 2
    for literals in c2_prime_term_literals:
        total_literals = prefix_literals + literals
        m3 += (total_literals - 1)

    # Part 4: a2b2'c2'
    m4 = 0
    prefix_literals = 2
    for literals in c2_prime_term_literals:
        total_literals = prefix_literals + literals
        m4 += (total_literals - 1)
        
    total_multiplications = m1 + m2 + m3 + m4

    print(f"Number of multiplications from a2'b2'c2 term: {m1}")
    print(f"Number of multiplications from a2b2c2 term: {m2}")
    print(f"Number of multiplications from a2'b2c2' term: {m3}")
    print(f"Number of multiplications from a2b2'c2' term: {m4}")
    print("-" * 20)
    print(f"Total multiplications = {m1} + {m2} + {m3} + {m4} = {total_multiplications}")

    print(f"\n<<<52>>>")

solve()