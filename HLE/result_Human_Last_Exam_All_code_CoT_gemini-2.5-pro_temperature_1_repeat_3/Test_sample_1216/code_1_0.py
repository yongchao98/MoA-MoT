import math

def solve_augustus_records():
    """
    Solves the puzzle about Augustus's imperial records based on the provided text and constraints.
    """
    # Step 1: Initial Counts
    total_records = 720
    lost_records = total_records // 6  # (a)
    readable_records = total_records - lost_records

    # Step 2: Use constraints and fractional clues to find the four key values.
    
    # (a) Lost records (primarily "Octavius")
    val_a = lost_records

    # (c) Single-variant documents ("Octavianus" or "Augustus")
    # Constraint: sqrt(s) = lost_records / 4
    # This leads to a number larger than the total readable records, but we follow the formula as written.
    val_c = (lost_records // 4) ** 2

    # (b) Dual-named documents ("Augustus" and "Caesar")
    # The constraint sqrt(d) implies 'd' must be a perfect square.
    # The most plausible interpretation of the fractional text is:
    # 1/3 of readable records are "Imperator Caesar Divi Filius" = 600/3 = 200
    # 1/4 of the *remaining* readable records use both names = (600 - 200) / 4 = 100.
    val_b = 100

    # (d) Full imperial title documents ("Imperator Caesar Augustus")
    # Constraint: f = sqrt(d) + 3 * l
    val_d = math.sqrt(val_b) + 3 * val_a

    # Step 3: Count distinct naming patterns for the divisor.
    # 1. Octavius (lost), 2. Imperator Caesar Divi Filius, 3. Augustus+Caesar,
    # 4. Octavianus/Augustus, 5. Imperator Caesar Augustus, 6. Caesar only
    distinct_patterns = 6

    # Step 4: (Optional but good for verification) Check divisibility constraint
    # Sum of all found numbers: l, d, s, f, plus others mentioned.
    # ICDF records = readable_records // 3 = 200
    # Caesar-only records = 80 (given)
    # sum_of_all = val_a + val_b + val_c + val_d + 200 + 80
    # print(f"Divisibility check: {sum_of_all} / {distinct_patterns} = {sum_of_all / distinct_patterns}")
    
    # Step 5: Calculate the final answer
    product = val_a * val_b * val_c * val_d
    division_result = product / distinct_patterns
    final_answer = math.floor(division_result / 1000)

    # Output the final equation with all its components
    print("Equation Breakdown:")
    print(f"  a) Lost records = {val_a}")
    print(f"  b) Dual-named documents = {int(val_b)}")
    print(f"  c) Single-variant documents = {val_c}")
    print(f"  d) Full imperial title documents = {int(val_d)}")
    print(f"  Divisor (Distinct naming patterns) = {distinct_patterns}")
    print("\nFinal Calculation:")
    print(f"  Product = {val_a} * {int(val_b)} * {val_c} * {int(val_d)} = {int(product)}")
    print(f"  Result = Product / Divisor = {int(product)} / {distinct_patterns} = {int(division_result)}")
    print(f"  Answer = floor(Result / 1000) = floor({int(division_result)} / 1000) = {final_answer}")
    
    # Return the final answer in the specified format
    print(f"\n<<< {final_answer} >>>")

solve_augustus_records()