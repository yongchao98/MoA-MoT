import sympy

def solve_triangle_ratio_problem():
    """
    This function provides a step-by-step derivation for the geometry problem
    and prints the final answer in terms of the triangle's side lengths.
    """

    # Define symbolic variables for the side lengths a, b, and c
    a, b, c = sympy.symbols('a b c')

    print("Derivation of the ratio BM/MI:")
    print("=" * 30)

    # Step 1: Explain the application of Ptolemy's Theorem
    print("1. The points A, B, C, and M all lie on the circumcircle, forming a cyclic quadrilateral ABCM.")
    print("2. According to Ptolemy's Theorem for a cyclic quadrilateral:")
    print("   (AB * CM) + (BC * AM) = (AC * BM)")
    print("   Substituting the side lengths a, b, c, we get:")
    print(f"   (c * CM) + (a * AM) = (b * BM)\n")

    # Step 2: Use the properties of point M
    print("3. The point M is where the angle bisector of angle B intersects the circumcircle.")
    print("   This means M bisects the arc AC, which implies the chords AM and CM are equal.")
    print("   Therefore, AM = CM.\n")

    # Step 3: Simplify the equation
    print("4. Substituting CM = AM into the Ptolemy's equation:")
    print("   (c * AM) + (a * AM) = (b * BM)")
    print("   AM * (a + c) = b * BM\n")

    # Step 4: Express the ratio of BM to AM
    print("5. Rearranging the equation gives the ratio of BM to AM:")
    print("   BM / AM = (a + c) / b\n")

    # Step 5: Apply the Trillium Theorem
    print("6. The Trillium Theorem (or Incenter-Excenter Lemma) states that for a point M defined this way,")
    print("   it is equidistant from A, C, and the incenter I. So, MA = MI = MC.")
    print("   This means we can replace AM with MI.\n")

    # Step 6: State the final result
    print("7. By substituting MI for AM, we arrive at the final ratio:")
    print("   BM / MI = (a + c) / b")
    
    print("=" * 30)
    print("Final Result:")
    
    # Define the final expression using sympy for pretty printing
    final_expression = (a + c) / b
    
    # Print the final equation
    print("The final expression for the ratio BM / MI is:")
    sympy.pprint(final_expression)
    
    print("\nIn the final equation, the components are:")
    print(f"Numerator: a + c (the sum of the sides adjacent to vertex B)")
    print(f"Denominator: b (the side opposite to vertex B)")

# Execute the function to display the solution
solve_triangle_ratio_problem()
