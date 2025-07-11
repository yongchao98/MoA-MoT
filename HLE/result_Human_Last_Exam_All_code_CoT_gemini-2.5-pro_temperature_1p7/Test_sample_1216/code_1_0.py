import math

def solve_augustus_records():
    """
    Solves the Augustus records puzzle based on the provided text and conditions.
    """
    # Step 1 & 2: Calculate 'a' (Lost records)
    total_records = 720
    a = total_records // 6

    # Step 3: Calculate 'c' (Single-variant documents) using the second condition
    # sqrt(c) = a / 4
    sqrt_c = a // 4
    c = sqrt_c ** 2

    # Step 4: Calculate 'd' (Full imperial title documents) using the ratio derived from the text
    # d/c = 3/4
    d = (3 * c) // 4

    # Step 5: Calculate 'b' (Dual-named documents) using the first condition
    # sqrt(b) + 3*a = d
    sqrt_b = d - 3 * a
    b = sqrt_b ** 2

    # Step 6: Define the product and the denominator for the final calculation
    product = a * b * c * d
    # The denominator is the sum of counts for patterns mentioning "Augustus"
    denominator = b + c + d

    # Step 7: Perform the final calculation
    result = product / denominator
    final_answer = math.floor(result / 1000)

    # Output the equation with the calculated numbers
    print("The final calculation is based on the equation:")
    print(f"(a * b * c * d) / (b + c + d) / 1000")
    print("Substituting the calculated values:")
    print(f"({a} * {b} * {c} * {d}) / ({b} + {c} + {d}) / 1000")
    print(f"Result: {final_answer}")

solve_augustus_records()