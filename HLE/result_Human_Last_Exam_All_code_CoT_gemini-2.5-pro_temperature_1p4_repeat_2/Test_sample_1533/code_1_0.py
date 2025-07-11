import sympy

def solve_geometry_problem():
    """
    This function solves the geometry problem by expressing the ratio BM/MI 
    in terms of the side lengths a, b, and c of triangle ABC using symbolic math.
    
    It requires the sympy library. You can install it using:
    pip install sympy
    """
    
    # --- Step 1: Define symbolic variables ---
    # a, b, c are the side lengths opposite to vertices A, B, C respectively.
    # s is the semi-perimeter of the triangle.
    a, b, c = sympy.symbols('a b c', positive=True)
    s = (a + b + c) / 2

    # --- Step 2: Explain the geometric reasoning ---
    print("Step 1: Deconstruct the ratio.")
    print("The points B, I, M are collinear, so BM = BI + IM.")
    print("The required ratio is BM / MI = (BI + IM) / IM = 1 + BI / IM.\n")

    print("Step 2: Use the Incenter-Circumcircle Lemma.")
    print("A key geometric property is that MA = MI, where M is the intersection of the angle bisector BI with the circumcircle.")
    print("Therefore, the ratio simplifies to: BM / MI = 1 + BI / MA.\n")
    
    print("Step 3: Express BI/MA using half-angle formulas in terms of side lengths.")
    print("The ratio BI / MA can be shown to be equal to (2 * sin(A/2) * sin(C/2)) / sin(B/2).")
    
    # Half-angle formulas for sine in terms of side lengths
    sin_A_half = sympy.sqrt((s - b) * (s - c) / (b * c))
    sin_B_half = sympy.sqrt((s - a) * (s - c) / (a * c))
    sin_C_half = sympy.sqrt((s - a) * (s - b) / (a * b))
    
    print("Using half-angle formulas:")
    print(f"  sin(A/2) = {sin_A_half}")
    print(f"  sin(B/2) = {sin_B_half}")
    print(f"  sin(C/2) = {sin_C_half}\n")

    # --- Step 3: Symbolically compute the ratio BI / MA ---
    print("Step 4: Compute the ratio BI / MA and simplify.")
    # The expression for the ratio in terms of angles
    ratio_BI_MA = (2 * sin_A_half * sin_C_half) / sin_B_half
    
    # Simplify the expression using sympy
    simplified_ratio_BI_MA = sympy.simplify(ratio_BI_MA)
    
    print(f"BI / MA = (2 * sin(A/2) * sin(C/2)) / sin(B/2) simplifies to:")
    print(f"BI / MA = {simplified_ratio_BI_MA}\n")
    
    # --- Step 4: Compute the final ratio BM / MI ---
    print("Step 5: Compute the final ratio BM / MI.")
    # The final ratio is 1 + BI/MA
    final_ratio = 1 + simplified_ratio_BI_MA
    
    # Simplify the final expression
    simplified_final_ratio = sympy.simplify(final_ratio)
    
    print("BM / MI = 1 + BI / MA")
    print(f"BM / MI = 1 + {simplified_ratio_BI_MA}")
    print(f"BM / MI = {simplified_final_ratio}\n")
    
    # --- Step 5: Print the final answer in the required format ---
    print("--- Final Answer ---")
    print("The expression for the ratio BM/MI in terms of side lengths a, b, and c is:")
    
    # As requested, outputting each "number" (symbol in this case) in the final equation
    numerator_vars = simplified_final_ratio.as_numer_denom()[0].args
    denominator_vars = simplified_final_ratio.as_numer_denom()[1]
    
    print(f"({numerator_vars[0]} + {numerator_vars[1]}) / {denominator_vars}")


if __name__ == '__main__':
    solve_geometry_problem()

<<<(a+c)/b>>>