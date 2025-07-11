import math

def solve_augustus_records():
    """
    Solves the puzzle about Augustus's imperial records.
    The solution is based on interpreting the proportional statements
    as a hierarchy, as a literal interpretation of all constraints leads to
    a logical contradiction. This path yields a consistent set of numbers
    that also satisfies the divisibility check mentioned in the prompt.
    """
    total_records = 720

    # Step 1: Calculate Lost Records (a) and Readable Records
    # "One-sixth of records were lost, primarily those using 'Octavius'"
    a = total_records // 6
    readable_records = total_records - a

    # Step 2: Separate the "Imperator Caesar Divi Filius" records
    # "One-third of remaining documents use 'Imperator Caesar Divi Filius'"
    divi_filius_records = readable_records // 3
    
    # Create the remaining pool for other categories
    pool = readable_records - divi_filius_records

    # Step 3: Sequentially calculate b, c, and d from the pool
    # "Among readable records: One-fourth use both 'Augustus' and 'Caesar'"
    # This is interpreted as 1/4 of the remaining pool.
    b = pool // 4

    # "Of those not using both names, 2/5 use either 'Octavianus' or 'Augustus'..."
    # "those not using both names" refers to the remainder of the pool after b is taken.
    pool_after_b = pool - b
    c = (pool_after_b * 2) // 5

    # "Half of remaining records use the full 'Imperator Caesar Augustus'"
    # "remaining records" refers to the remainder of the pool after b and c are taken.
    pool_after_c = pool_after_b - c
    d = pool_after_c // 2
    
    print(f"Lost records (a): {a}")
    print(f"Dual-named documents (b): {b}")
    print(f"Single-variant documents (c): {c}")
    print(f"Full imperial title documents (d): {d}")
    print("-" * 20)

    # Step 4: Calculate the product
    product = a * b * c * d
    
    # Step 5: Define the divisor
    # The distinct patterns are: Divi Filius, dual-named(b), single-variant(c),
    # full title(d), Caesar only(e=80), and the remainder(g). This is 6 patterns.
    # The sum of all found numbers (720) is divisible by 6, confirming this.
    divisor = 6
    
    print("Final Equation Calculation:")
    print(f"Product = a * b * c * d = {a} * {b} * {c} * {d} = {product}")
    print(f"Divisor (Number of distinct readable patterns) = {divisor}")

    # Step 6: Calculate the result and the final answer
    result = product / divisor
    final_answer = math.floor(result / 1000)

    print(f"Result = Product / Divisor = {product} / {divisor} = {int(result)}")
    print(f"Final Answer = floor(Result / 1000) = floor({int(result)} / 1000) = {final_answer}")
    
    # The final output required by the prompt
    print(f"\n<<<{final_answer}>>>")

solve_augustus_records()