import math

def solve_augustus_records_puzzle():
    """
    This function solves the word puzzle about Augustus' imperial records
    by calculating the values of a, b, c, and d based on the given constraints.
    """

    # Step 1: Define initial values from the problem statement.
    total_records = 720
    # The number of records using only "Caesar" is given.
    caesar_only_records = 80
    # The number of "Imperator Caesar Divi Filius" records.
    # One-third of the remaining (readable) records.
    readable_records = total_records * (5/6)
    icdf_records = math.floor(readable_records / 3)

    # Step 2: Calculate 'a' (Lost records, primarily "Octavius").
    # "One-sixth of records were lost"
    a = math.floor(total_records / 6)

    # Step 3: Calculate 'c' (Single-variant documents).
    # "The records using single variants...form a perfect square whose root equals
    # the number of lost records divided by 4"
    # This reveals the core puzzle logic, as c > readable_records.
    sqrt_c = a / 4
    c = int(sqrt_c ** 2)

    # Step 4: Calculate 'd' (Full imperial title documents).
    # This requires interpreting the chained proportional statements.
    # d = 1/2 * ( (K-b) - c )
    # c = 2/5 * (K-b)  =>  (K-b) = 5/2 * c
    # By substituting (K-b), we can solve for d without knowing K or b.
    k_minus_b = (5/2) * c
    # "Half of remaining records use the full 'Imperator Caesar Augustus'"
    d = math.floor(0.5 * (k_minus_b - c))

    # Step 5: Calculate 'b' (Dual-named "Augustus" and "Caesar" records).
    # "The square root of dual-named records...plus thrice the Octavius-lost records
    # equals the number of records using the full imperial title"
    # sqrt(b) + 3*a = d  => sqrt(b) = d - 3*a
    sqrt_b = d - 3 * a
    b = int(sqrt_b ** 2)

    # Step 6: Define the distinct naming patterns and their counts for the denominator.
    # The denominator is the sum of the counts for each distinct way Augustus is named.
    patterns = {
        "Octavius (Lost)": a,
        "Imperator Caesar Divi Filius": icdf_records,
        "Augustus and Caesar (Dual-named)": b,
        "Octavianus or Augustus (Single-variant)": c,
        "Imperator Caesar Augustus (Full Title)": d,
        "Only Caesar": caesar_only_records
    }
    denominator = sum(patterns.values())

    # Step 7: Calculate the numerator (product of a, b, c, d).
    numerator = a * b * c * d
    
    # Step 8: Perform the final calculation as requested.
    answer = numerator / denominator
    final_answer = math.floor(answer / 1000)

    # Step 9: Print the results clearly.
    print("--- Calculated Values ---")
    print(f"a) Lost records (Octavius ones): {a}")
    print(f"b) Dual-named (Augustus+Caesar) documents: {b}")
    print(f"c) Single-variant (Octavianus/Augustus) documents: {c}")
    print(f"d) Full imperial title documents: {d}")
    print("\n--- Final Calculation ---")
    print("The final equation is: floor( (a * b * c * d) / (Sum of records for all distinct patterns) / 1000 )")
    print(f"The equation with numbers is: floor( ({a} * {b} * {c} * {d}) / {denominator} / 1000 )")
    print(f"Result of the division: {answer:.2f}")

    print("\nFinal Answer:")
    print(final_answer)
    
solve_augustus_records_puzzle()
<<<71559>>>