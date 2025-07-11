import math

def calculate_l(d: int):
    """
    Calculates the exact value of l(d) based on the derived formula.

    The formula is l(d) = 1 + ln((sqrt(d)-1)/(sqrt(d)+1)).

    Args:
        d: An integer greater than or equal to 2.
    """
    if not isinstance(d, int) or d < 2:
        print("Error: d must be an integer greater than or equal to 2.")
        return

    print(f"Calculating l(d) for d = {d}")
    
    # Step 1: Calculate sqrt(d)
    sqrt_d = math.sqrt(d)
    print(f"sqrt(d) = {sqrt_d}")

    # Step 2: Calculate the numerator of the ratio
    numerator = sqrt_d - 1
    print(f"sqrt(d) - 1 = {numerator}")

    # Step 3: Calculate the denominator of the ratio
    denominator = sqrt_d + 1
    print(f"sqrt(d) + 1 = {denominator}")

    # Step 4: Calculate the ratio
    ratio = numerator / denominator
    print(f"(sqrt(d) - 1) / (sqrt(d) + 1) = {ratio}")

    # Step 5: Calculate the natural logarithm of the ratio
    log_val = math.log(ratio)
    print(f"ln(ratio) = {log_val}")

    # Step 6: Calculate the final result
    result = 1 + log_val
    
    # Final equation output
    print("\nFinal Equation:")
    print(f"l({d}) = 1 + ln((sqrt({d}) - 1) / (sqrt({d}) + 1))")
    print(f"l({d}) = 1 + ln(({numerator}) / ({denominator}))")
    print(f"l({d}) = 1 + ln({ratio})")
    print(f"l({d}) = 1 + ({log_val})")
    print(f"l({d}) = {result}")


# Example usage with d=4
calculate_l(4)

# Example usage with d=9
# calculate_l(9)