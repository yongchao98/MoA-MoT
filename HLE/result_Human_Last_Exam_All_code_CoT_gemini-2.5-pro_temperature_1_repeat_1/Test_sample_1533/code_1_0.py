import sympy

def solve_geometry_ratio():
    """
    Solves the triangle geometry problem symbolically using sympy.
    """
    # Define the symbols for the side lengths of triangle ABC and other lengths.
    a, b, c, BM, MI, MA, MC = sympy.symbols('a b c BM MI MA MC')

    # Step 1: State the key geometric properties as equations.
    # Property 1: M is equidistant from A and C because it's on the angle bisector of B.
    # Property 2: M is equidistant from A and I (Incenter-Excenter Lemma).
    # This gives us MI = MA = MC.
    eq_ma_mc = sympy.Eq(MA, MC)
    eq_ma_mi = sympy.Eq(MA, MI)
    print("Step 1: Establish key equalities from the geometry of the problem.")
    print(f"From the properties of point M on the circumcircle, we know that MA = MC and MA = MI.")
    print(f"Represented symbolically: {eq_ma_mc} and {eq_ma_mi}\n")

    # Step 2: Apply Ptolemy's theorem to the cyclic quadrilateral ABCM.
    # The theorem states: AB * CM + BC * AM = AC * BM
    print("Step 2: Apply Ptolemy's Theorem to the cyclic quadrilateral ABCM.")
    print("The theorem states: AB * CM + BC * AM = AC * BM")
    eq_ptolemy = sympy.Eq(c * MC + a * MA, b * BM)
    print(f"In terms of a, b, c: {eq_ptolemy}\n")

    # Step 3: Substitute the equalities from Step 1 into Ptolemy's equation.
    # Substitute MC with MA.
    print("Step 3: Substitute the equalities into Ptolemy's equation.")
    eq_substituted = eq_ptolemy.subs(MC, MA)
    print(f"Substituting MC with MA gives: {eq_substituted}")

    # Factor out MA
    eq_factored = sympy.Eq(sympy.factor(eq_substituted.lhs), eq_substituted.rhs)
    print(f"Factoring out MA gives: {eq_factored}\n")

    # Step 4: Solve for the ratio BM/MA.
    print("Step 4: Solve for the ratio BM / MA.")
    # To solve for the ratio, we can express BM in terms of MA.
    bm_expr = sympy.solve(eq_factored, BM)[0]
    # The ratio is BM / MA
    ratio_bm_ma = bm_expr / MA
    print(f"The ratio BM / MA is: {ratio_bm_ma}\n")

    # Step 5: Substitute MA with MI to get the final desired ratio BM/MI.
    print("Step 5: Since MA = MI, the final ratio BM / MI is the same.")
    final_ratio = ratio_bm_ma
    
    # The instruction asks to output each number in the final equation.
    # We can represent the final expression as a string to highlight the coefficients.
    # For (a+c)/b, the coefficients are 1.
    final_expression_str = f"(1*a + 1*c) / (1*b)"

    print("-----------------------------------------")
    print("Final Answer")
    print("-----------------------------------------")
    print(f"The expression for the ratio BM / MI is:")
    sympy.pprint(final_ratio)
    #print(f"\nTo explicitly show the numbers in the equation as requested:")
    #print(final_expression_str)


if __name__ == '__main__':
    solve_geometry_ratio()