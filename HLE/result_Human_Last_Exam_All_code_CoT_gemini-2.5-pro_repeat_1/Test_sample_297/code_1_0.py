def count_multiplications():
    """
    Calculates the number of multiplication operations in the fully expanded
    expression for the sum bit s2 of a 3-bit binary adder.
    """
    
    # 1. Define the minimal SOP for the carry bit c2.
    # c2 = a1b1 + a0a1b0 + a0b1b0
    # We represent each product term by the number of its literals.
    c2_literal_counts = [2, 4, 4]
    
    # 2. Define the minimal SOP for the complement of the carry bit, c2'.
    # c2' = a1'b1' + a1'a0' + a1'b0' + b1'a0' + b1'b0'
    # Each term has 2 literals.
    c2_prime_literal_counts = [2, 2, 2, 2, 2]
    
    total_multiplications = 0
    
    # 3. The SOP expression for s2 is: a2'b2'c2 + a2'b2c2' + a2b2'c2' + a2b2c2
    # We calculate the multiplications for each of these four components.
    
    # Component 1: a2'b2'c2
    # Expanding this means multiplying a2'b2' with each term of c2.
    # This adds 2 literals to each term in c2_literal_counts.
    # A term with 'n' literals has 'n-1' multiplications.
    comp1_mults = 0
    for literals in c2_literal_counts:
        final_literals = literals + 2  # Add literals from a2'b2'
        comp1_mults += (final_literals - 1)
    
    # Component 2: a2'b2c2'
    # This adds 2 literals (a2', b2) to each term in c2_prime_literal_counts.
    comp2_mults = 0
    for literals in c2_prime_literal_counts:
        final_literals = literals + 2
        comp2_mults += (final_literals - 1)
        
    # Component 3: a2b2'c2'
    # This adds 2 literals (a2, b2') to each term in c2_prime_literal_counts.
    comp3_mults = 0
    for literals in c2_prime_literal_counts:
        final_literals = literals + 2
        comp3_mults += (final_literals - 1)
        
    # Component 4: a2b2c2
    # This adds 2 literals (a2, b2) to each term in c2_literal_counts.
    comp4_mults = 0
    for literals in c2_literal_counts:
        final_literals = literals + 2
        comp4_mults += (final_literals - 1)
        
    # 4. Sum the multiplications from all components for the final answer.
    total_multiplications = comp1_mults + comp2_mults + comp3_mults + comp4_mults
    
    print("The total number of multiplication operations is derived by summing the operations from each part of the s2 expression:")
    print(f"Multiplications from a2'b2'c2 = {comp1_mults}")
    print(f"Multiplications from a2'b2c2' = {comp2_mults}")
    print(f"Multiplications from a2b2'c2' = {comp3_mults}")
    print(f"Multiplications from a2b2c2  = {comp4_mults}")
    print(f"Total multiplications = {comp1_mults} + {comp2_mults} + {comp3_mults} + {comp4_mults} = {total_multiplications}")

count_multiplications()