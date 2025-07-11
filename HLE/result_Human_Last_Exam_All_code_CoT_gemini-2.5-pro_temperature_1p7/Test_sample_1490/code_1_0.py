from fractions import Fraction

def solve_rangoli_riddle():
    """
    Calculates the total number of curves in the restored rangoli pattern
    based on the poem's scenario.
    """
    # 1. Start with the known value from the poem.
    curves_stayed_true = 90
    print(f"The number of curves that remained unchanged is {curves_stayed_true}.")
    print("-" * 20)

    # 2. Define the fractions of disturbed curves from the poem.
    frac_lost_shape = Fraction(3, 8)
    frac_new_path = Fraction(1, 4)

    # 3. Calculate the fraction of curves that stayed true.
    # The total (1) minus the disturbed fractions gives the fraction that remained.
    frac_stayed_true = 1 - frac_lost_shape - frac_new_path
    
    print(f"Fraction of curves that lost shape: {frac_lost_shape}")
    print(f"Fraction of curves that found new paths: {frac_new_path}")
    print(f"Therefore, the fraction of curves that stayed true is 1 - {frac_lost_shape} - {frac_new_path} = {frac_stayed_true}.")
    print("-" * 20)
    
    # 4. Calculate the original total number of curves.
    # If `frac_stayed_true` of the total is 90, we can find the total.
    total_curves = int(curves_stayed_true / frac_stayed_true)
    print(f"The original total number of curves is calculated as: {curves_stayed_true} / {frac_stayed_true} = {total_curves}.")
    print("-" * 20)

    # 5. Calculate the number of curves in each category.
    curves_lost_shape = int(total_curves * frac_lost_shape)
    curves_new_path = int(total_curves * frac_new_path)
    
    # Optional: Calculate the breakdown of the transformed curves to show full understanding.
    parabolic = int(curves_lost_shape * Fraction(1, 5))
    elliptical = int(curves_lost_shape * Fraction(2, 9))
    circular = curves_lost_shape - parabolic - elliptical

    print("Breaking down the final pattern:")
    print(f"- Curves that stayed true: {curves_stayed_true}")
    print(f"- Curves that were transformed (lost shape): {curves_lost_shape}")
    print(f"    - Parabolic: {parabolic}")
    print(f"    - Elliptical: {elliptical}")
    print(f"    - Circular: {circular}")
    print(f"- Curves that found new pathways: {curves_new_path}")
    print("-" * 20)

    # 6. The total number of curves in the restored pattern is the sum of its parts.
    # The golden ratio affects the properties, not the count, of the new curves.
    print("To restore the pattern, the master must account for all curves.")
    print("The final total is the sum of the unchanged, transformed, and new pathway curves:")
    print(f"{curves_stayed_true} + {curves_lost_shape} + {curves_new_path} = {total_curves}")

solve_rangoli_riddle()
<<<240>>>