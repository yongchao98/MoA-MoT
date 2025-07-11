def calculate_multiplications():
    """
    Calculates and explains the number of multiplications in the
    fully expanded expression for the sum bit s2 of a 3-bit binary adder.
    """
    print("This script calculates the number of multiplication operations in the")
    print("fully expanded boolean expression for s2, the sum bit at position 2,")
    print("when adding two 3-bit binary numbers A=a2a1a0 and B=b2b1b0.")
    print("\n---")

    # Explain the formula for s2
    print("The sum bit s2 is given by the formula: s2 = a2 ⊕ b2 ⊕ c2,")
    print("where c2 is the carry-in from the previous bit (bit 1).")
    print("In the requested sum-of-products form (using only addition, multiplication, and single-digit complement), this is expressed as:")
    print("s2 = (a2'b2'c2) + (a2b2c2) + (a2'b2c2') + (a2b2'c2')")
    print("\n---")

    # Explain c2 and c2'
    print("The carry c2 is generated from adding a1, b1, and the carry c1 (where c1 = a0*b0).")
    print("The fully expanded expression for c2 is: c2 = a1*b1 + a1*a0*b0 + b1*a0*b0")
    print("The fully expanded expression for its complement c2' is: c2' = a1'*b1' + a1'*a0' + a1'*b0' + b1'*a0' + b1'*b0'")
    print("\n---")

    print("We now count the multiplications in the full expansion of s2 part by part.")
    print("A product of 'n' literals (e.g., a*b*c) requires n-1 multiplications.")
    print("\n---")

    # Part 1: Expansion of a2'b2'c2
    print("Part 1: a2'*b2'*c2 = a2'*b2'*(a1*b1 + a1*a0*b0 + b1*a0*b0)")
    print("This expands to three product terms:")
    mults_1_1 = 4 - 1
    print(f"  - Term a2'*b2'*a1*b1 has 4 literals, requiring {mults_1_1} multiplications.")
    mults_1_2 = 5 - 1
    print(f"  - Term a2'*b2'*a1*a0*b0 has 5 literals, requiring {mults_1_2} multiplications.")
    mults_1_3 = 5 - 1
    print(f"  - Term a2'*b2'*b1*a0*b0 has 5 literals, requiring {mults_1_3} multiplications.")
    subtotal_1 = mults_1_1 + mults_1_2 + mults_1_3
    print(f"The equation for the subtotal of Part 1 is: {mults_1_1} + {mults_1_2} + {mults_1_3} = {subtotal_1}")
    print("---")

    # Part 2: Expansion of a2b2c2
    print("Part 2: a2*b2*c2 = a2*b2*(a1*b1 + a1*a0*b0 + b1*a0*b0)")
    print("This part is symmetric to Part 1 and also expands to three terms:")
    mults_2_1 = 4 - 1
    print(f"  - Term a2*b2*a1*b1 has 4 literals, requiring {mults_2_1} multiplications.")
    mults_2_2 = 5 - 1
    print(f"  - Term a2*b2*a1*a0*b0 has 5 literals, requiring {mults_2_2} multiplications.")
    mults_2_3 = 5 - 1
    print(f"  - Term a2*b2*b1*a0*b0 has 5 literals, requiring {mults_2_3} multiplications.")
    subtotal_2 = mults_2_1 + mults_2_2 + mults_2_3
    print(f"The equation for the subtotal of Part 2 is: {mults_2_1} + {mults_2_2} + {mults_2_3} = {subtotal_2}")
    print("---")

    # Part 3: Expansion of a2'b2c2'
    print("Part 3: a2'*b2*c2' = a2'*b2*(a1'*b1' + a1'*a0' + a1'*b0' + b1'*a0' + b1'*b0')")
    num_terms_3 = 5
    literals_per_term_3 = 4 # 2 literals from a2'*b2 and 2 from each term in c2'
    mults_per_term_3 = literals_per_term_3 - 1
    subtotal_3 = num_terms_3 * mults_per_term_3
    print(f"This expands to {num_terms_3} terms, each with {literals_per_term_3} literals.")
    print(f"Each term requires {mults_per_term_3} multiplications.")
    print(f"The equation for the subtotal of Part 3 is: {num_terms_3} * {mults_per_term_3} = {subtotal_3}")
    print("---")
    
    # Part 4: Expansion of a2b2'c2'
    print("Part 4: a2*b2'*c2' = a2*b2'*(a1'*b1' + a1'*a0' + a1'*b0' + b1'*a0' + b1'*b0')")
    num_terms_4 = 5
    literals_per_term_4 = 4 # 2 literals from a2*b2' and 2 from each term in c2'
    mults_per_term_4 = literals_per_term_4 - 1
    subtotal_4 = num_terms_4 * mults_per_term_4
    print(f"This part is symmetric to Part 3 and results in {num_terms_4} terms with {literals_per_term_4} literals each.")
    print(f"Each term requires {mults_per_term_4} multiplications.")
    print(f"The equation for the subtotal of Part 4 is: {num_terms_4} * {mults_per_term_4} = {subtotal_4}")
    print("---")
    
    # Final total
    total_mults = subtotal_1 + subtotal_2 + subtotal_3 + subtotal_4
    print("\nFinal Calculation:")
    print("The final equation for the total number of multiplications is the sum of the subtotals:")
    print(f"{subtotal_1} + {subtotal_2} + {subtotal_3} + {subtotal_4} = {total_mults}")
    
    print(f"\nThus, the total number of multiplication operations in the fully expanded expression for s2 is {total_mults}.")

calculate_multiplications()