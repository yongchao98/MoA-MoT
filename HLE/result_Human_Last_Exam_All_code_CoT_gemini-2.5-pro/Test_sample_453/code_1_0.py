import fractions

def calculate_relativistic_correction():
    """
    Calculates the relativistic energy correction for a hydrogen atom state.
    
    The function substitutes the given quantum numbers into the standard formula
    for the first-order relativistic energy correction and prints the step-by-step
    calculation.
    """
    # Given quantum numbers
    n = 3
    l = 2

    print(f"Calculating the relativistic energy shift for n = {n}, l = {l}.")
    print("-" * 30)
    
    # The formula for the relativistic energy correction is:
    # ΔE = - (mc²α⁴ / (2n⁴)) * (n/(l+1/2) - 3/4)
    print("The formula for the energy shift is:")
    print("ΔE = - (m⋅c²⋅α⁴ / (2⋅n⁴)) * (n/(l + 1/2) - 3/4)")
    print("\nSubstituting the values n=3 and l=2:")
    
    # Step 1: Substitute n and l into the formula string
    print(f"ΔE = - (m⋅c²⋅α⁴ / (2⋅{n}⁴)) * ({n}/({l} + 1/2) - 3/4)")
    
    # Step 2: Evaluate the denominators
    print(f"ΔE = - (m⋅c²⋅α⁴ / (2⋅{n**4})) * ({n}/(2.5) - 3/4)")
    
    # Step 3: Continue simplifying
    print(f"ΔE = - (m⋅c²⋅α⁴ / {2*n**4}) * (6/5 - 3/4)")
    
    # Step 4: Find a common denominator for the terms in the parenthesis
    print(f"ΔE = - (m⋅c²⋅α⁴ / 162) * ((24 - 15)/20)")
    
    # Step 5: Simplify the parenthesis
    print(f"ΔE = - (m⋅c²⋅α⁴ / 162) * (9/20)")
    
    # Step 6: Combine the fractions
    print(f"ΔE = - (9 ⋅ m⋅c²⋅α⁴) / (162 ⋅ 20)")
    
    # Step 7: Simplify the final fraction
    # 162 / 9 = 18
    print(f"ΔE = - (m⋅c²⋅α⁴) / (18 ⋅ 20)")
    
    # Final result
    final_denominator = (162 * 20) // 9
    print(f"\nFinal Result:")
    print(f"ΔE = - (m⋅c²⋅α⁴) / {final_denominator}")

# Run the calculation
calculate_relativistic_correction()