import sympy

def solve_ratio_calculation():
    """
    This function demonstrates the calculation of the ratio BM/MI for a sample triangle.

    The derivation is as follows:
    1.  Let M be the intersection of the angle bisector of B with the circumcircle. By the Trillium Theorem (or Incenter-Excenter Lemma), we can prove that MA = MC = MI.
        - Proof of MA = MC: The angle bisector BM splits the arc AC into two equal arcs, AM and MC. Equal arcs mean equal chords, so MA = MC.
        - Proof of MA = MI: We show that triangle MAI is isosceles.
          - Angle MIA is the exterior angle to triangle AIB, so ∠MIA = ∠IAB + ∠IBA = A/2 + B/2.
          - Angle MAI = ∠MAC + ∠CAI. Since M is on the circumcircle, ∠MAC = ∠MBC = B/2. And AI is the angle bisector, so ∠CAI = A/2. Thus, ∠MAI = B/2 + A/2.
          - Since ∠MIA = ∠MAI, triangle MAI is isosceles with MA = MI.
    2.  The points A, B, C, M form a cyclic quadrilateral ABCM. By Ptolemy's Theorem:
        AB * CM + BC * AM = AC * BM
    3.  Substituting side lengths a, b, c and using AM = CM = MI:
        c * MI + a * MI = b * BM
        MI * (a + c) = b * BM
    4.  Solving for the ratio BM/MI:
        BM / MI = (a + c) / b
    """

    # We use symbolic variables to represent the sides
    a, b, c = sympy.symbols('a b c')

    # The derived ratio is (a+c)/b
    ratio_expr = (a + c) / b

    print("The derived ratio BM / MI is expressed in terms of side lengths a, b, and c.")
    print("Let a = BC, b = AC, c = AB.")
    print("The final derived formula is:")
    # The print function will output the final symbolic expression
    final_equation_str = f"BM / MI = (a + c) / b"
    print(final_equation_str)

    # Example calculation for a triangle with sides a=7, b=8, c=9
    a_val, b_val, c_val = 7, 8, 9
    numerator = a_val + c_val
    denominator = b_val
    result = numerator / denominator

    print("\n--- Example Calculation ---")
    print(f"For a triangle with side lengths a = {a_val}, b = {b_val}, c = {c_val}:")
    # In the final code you still need to output each number in the final equation!
    calculation_str = f"BM / MI = ({a_val} + {c_val}) / {b_val} = {numerator} / {denominator} = {result}"
    print(calculation_str)

# Run the function
solve_ratio_calculation()