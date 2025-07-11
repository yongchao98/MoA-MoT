import math

def solve_multiplication_count():
    """
    This script calculates the number of multiplications in the fully expanded
    boolean expression for the s2 sum bit of a 3-bit binary adder.
    
    The final expression is derived from s2 = a2 XOR b2 XOR c2, fully expanded
    into a sum-of-products form using only the input digits (a2,a1,a0, b2,b1,b0).
    """

    # 1. Define the sum-of-products terms for the carry bit c2.
    # c2 = a1*b1 + a1*a0*b0 + b1*a0*b0
    c2_terms_literals = [
        ['a1', 'b1'],         # 2 literals
        ['a1', 'a0', 'b0'],   # 3 literals
        ['b1', 'a0', 'b0']    # 3 literals
    ]

    # 2. Define the sum-of-products terms for the complement of the carry bit, c2'.
    # c2' = a1'*b1' + a1'*a0' + a1'*b0' + b1'*a0' + b1'*b0'
    c2_prime_terms_literals = [
        ["a1'", "b1'"],      # 2 literals
        ["a1'", "a0'"],      # 2 literals
        ["a1'", "b0'"],      # 2 literals
        ["b1'", "a0'"],      # 2 literals
        ["b1'", "b0'"]       # 2 literals
    ]

    total_multiplications = 0
    final_equation_terms = []

    # The expression for s2 is a2'b2'c2 + a2'b2c2' + a2b2'c2' + a2b2c2.
    # We will build each of these four parts.

    # Part 1: a2'b2'c2
    prefix_literals = ["a2'", "b2'"]
    for term in c2_terms_literals:
        # Combine prefix literals with c2 term literals
        full_term = prefix_literals + term
        final_equation_terms.append(" * ".join(full_term))
        # N literals require N-1 multiplications
        total_multiplications += len(full_term) - 1

    # Part 2: a2'b2c2'
    prefix_literals = ["a2'", "b2"]
    for term in c2_prime_terms_literals:
        full_term = prefix_literals + term
        final_equation_terms.append(" * ".join(full_term))
        total_multiplications += len(full_term) - 1
        
    # Part 3: a2b2'c2'
    prefix_literals = ["a2", "b2'"]
    for term in c2_prime_terms_literals:
        full_term = prefix_literals + term
        final_equation_terms.append(" * ".join(full_term))
        total_multiplications += len(full_term) - 1

    # Part 4: a2b2c2
    prefix_literals = ["a2", "b2"]
    for term in c2_terms_literals:
        full_term = prefix_literals + term
        final_equation_terms.append(" * ".join(full_term))
        total_multiplications += len(full_term) - 1
    
    print("The fully expanded expression for s2 is:\n")
    # Print the full equation by joining all the product terms with " + "
    final_equation_str = "s2 = " + " + \\\n     ".join(final_equation_terms)
    print(final_equation_str)
    
    print("\n--------------------------------------------------")
    print("Calculation of Total Multiplication Operations")
    print("--------------------------------------------------")
    # Calculating parts again just for explanation.
    # Part 1 (a2'b2'c2): 3 terms. (2+2-1) + (2+3-1) + (2+3-1) = 3 + 4 + 4 = 11
    part1_mults = (2+2-1) + (2+3-1) + (2+3-1)
    # Part 2 (a2'b2c2'): 5 terms. (2+2-1) * 5 = 3 * 5 = 15
    part2_mults = (2+2-1) * 5
    # Part 3 (a2b2'c2'): 5 terms. (2+2-1) * 5 = 3 * 5 = 15
    part3_mults = (2+2-1) * 5
    # Part 4 (a2b2c2): 3 terms. (2+2-1) + (2+3-1) + (2+3-1) = 3 + 4 + 4 = 11
    part4_mults = (2+2-1) + (2+3-1) + (2+3-1)
    
    print(f"Term 1 (a2'b2'c2) contributes: {part1_mults} multiplications")
    print(f"Term 2 (a2'b2c2') contributes: {part2_mults} multiplications")
    print(f"Term 3 (a2b2'c2') contributes: {part3_mults} multiplications")
    print(f"Term 4 (a2b2c2) contributes: {part4_mults} multiplications")
    print("--------------------------------------------------")

    print(f"The total number of multiplication operations is {total_multiplications}.")

solve_multiplication_count()
<<<52>>>