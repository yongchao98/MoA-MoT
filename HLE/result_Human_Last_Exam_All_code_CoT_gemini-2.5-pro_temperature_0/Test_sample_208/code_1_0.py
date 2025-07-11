import math

def calculate_packing_radius():
    """
    Calculates the radius of a circle that tightly packs fourteen circles of radius one.
    """
    # The radius of the small circles
    r = 1.0

    # The formula for the radius R of the large circle for N=14 small circles is:
    # R/r = 1 + sqrt(2 + 4/sqrt(3))
    # We will calculate this step by step.

    print(f"The radius R is found by the equation: R = r * (1 + sqrt(2 + 4 / sqrt(3)))")
    print(f"Given the small circle radius r = {r}\n")

    # Step 1: Calculate the value of sqrt(3)
    val_sqrt3 = math.sqrt(3)
    print(f"Step 1: Calculate sqrt(3)")
    print(f"sqrt(3) = {val_sqrt3}\n")

    # Step 2: Calculate the term 4 / sqrt(3)
    term1 = 4 / val_sqrt3
    print(f"Step 2: Calculate the term 4 / sqrt(3)")
    print(f"4 / {val_sqrt3} = {term1}\n")

    # Step 3: Calculate the term inside the main square root
    term2 = 2 + term1
    print(f"Step 3: Calculate the term inside the square root: 2 + (the result from Step 2)")
    print(f"2 + {term1} = {term2}\n")

    # Step 4: Take the square root of the result from Step 3
    term3 = math.sqrt(term2)
    print(f"Step 4: Take the square root of the result from Step 3")
    print(f"sqrt({term2}) = {term3}\n")

    # Step 5: Complete the calculation for R
    R = r * (1 + term3)
    print(f"Step 5: Complete the final calculation for R = r * (1 + result from Step 4)")
    print(f"R = {r} * (1 + {term3}) = {R}\n")

    # Round the result to 4 significant digits
    # We use the .4g format specifier for 4 significant digits
    rounded_R = f"{R:.4g}"

    print("---")
    print(f"The final calculated radius is: {R}")
    print(f"The radius of the circle up to 4 significant digits is: {rounded_R}")

calculate_packing_radius()