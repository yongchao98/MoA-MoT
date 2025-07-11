import sympy as sp

def solve_geometry_ratio():
    """
    This function symbolically derives the ratio BM/MI in terms of the side
    lengths a, b, and c of triangle ABC.
    """
    # Define symbolic variables for the side lengths
    a, b, c = sp.symbols('a b c', positive=True)

    # Define the semi-perimeter s
    s = (a + b + c) / 2

    # The ratio BM/MI can be expressed as BI/MI + 1.
    # We derived that BI/MI simplifies to 2*(s-b)/b.
    # Let's define this term.
    ratio_bi_mi = 2 * (s - b) / b
    
    print("Step 1: The ratio BM/MI is equal to BI/MI + 1.")
    print("Step 2: We derived the expression for BI/MI in terms of the semi-perimeter 's' and side 'b'.")
    print(f"BI / MI = {ratio_bi_mi}")
    
    # Simplify the expression for BI/MI
    simplified_bi_mi = sp.simplify(ratio_bi_mi)
    print("\nStep 3: Simplifying the expression for BI/MI in terms of a, b, c:")
    print(f"BI / MI = {simplified_bi_mi}")

    # Now, calculate the final ratio BM/MI = BI/MI + 1
    final_ratio = simplified_bi_mi + 1
    print("\nStep 4: Calculating the final ratio BM/MI = (BI/MI) + 1:")
    print(f"BM / MI = {simplified_bi_mi} + 1")

    # Simplify the final expression
    simplified_final_ratio = sp.simplify(final_ratio)
    
    # To satisfy the prompt "output each number in the final equation",
    # we will print the numerator and denominator of the final fraction.
    num, den = simplified_final_ratio.as_numer_denom()

    print("\n--- Final Result ---")
    print("The simplified expression for the ratio BM/MI is:")
    print(f"Numerator: {num}")
    print(f"Denominator: {den}")
    print(f"BM / MI = {simplified_final_ratio}")

if __name__ == '__main__':
    solve_geometry_ratio()
