import math

def solve_augustus_records():
    """
    Solves the problem by setting up and solving a system of equations
    derived from the problem's constraints.
    """
    total_records = 720
    lost_records_total = total_records // 6
    readable_records = total_records - lost_records_total

    # From the text, we can determine the counts of two distinct partitions of readable records:
    # f: "Imperator Caesar Divi Filius" records
    # e: "Caesar" only records
    f = readable_records // 3
    e = 80

    # The sum of all distinct partitions of readable records is 600.
    # The main partitions are f, e, c (single-variant), and b (dual-named).
    # This gives the equation: f + e + b + c = 600
    # 200 + 80 + b + c = 600  => b + c = 320
    b_plus_c = 320
    
    a = 0
    b = 0
    # We need to find 'a' (Octavius-lost records) that satisfies the constraints.
    # Constraint: c = (a/4)^2, which means c = a^2 / 16. 'a' must be a multiple of 4.
    # From b + c = 320, we get b = 320 - c = 320 - (a^2 / 16).
    # Constraint: 'b' must be a perfect square.
    # We iterate through multiples of 4 for 'a' to find a value that makes 'b' a perfect square.
    for val_a in range(4, 100, 4):
        val_b_float = b_plus_c - (val_a**2 / 16)
        if val_b_float > 0:
            val_b_sqrt = math.sqrt(val_b_float)
            if val_b_sqrt == int(val_b_sqrt):
                a = val_a
                b = int(val_b_float)
                break
    
    # Now that we have a and b, we can find c and d.
    # c = b + c - b
    c = b_plus_c - b

    # d is found from the constraint: d = sqrt(b) + 3*a
    d = int(math.sqrt(b)) + 3 * a

    # The problem asks for the product of a, b, c, d
    product = a * b * c * d

    # The denominator is "The sum of the distinct ways Augustus is named in the archive".
    # These are the records in partitions b and c.
    denominator = b + c

    # Calculate the result
    result = product / denominator

    # The final answer is floor(result / 1000)
    final_answer = math.floor(result / 1000)

    # Output the required values and the final equation.
    print(f"Based on the constraints, the calculated values are:")
    print(f"a) Lost records (Octavius ones): {a}")
    print(f"b) Documents with dual naming (Augustus+Caesar): {b}")
    print(f"c) Single-variant documents (Octavianus/Augustus): {c}")
    print(f"d) Full imperial title documents: {d}")
    print("\nCalculating the final answer:")
    print(f"Product = {a} * {b} * {c} * {d} = {product}")
    print(f"Denominator (sum of ways Augustus is named) = b + c = {b} + {c} = {denominator}")
    print(f"Result = Product / Denominator = {product} / {denominator} = {result}")
    print(f"Final Answer = floor(Result / 1000) = floor({result} / 1000) = {final_answer}")
    
    print(f"\n<<< {final_answer} >>>")

solve_augustus_records()