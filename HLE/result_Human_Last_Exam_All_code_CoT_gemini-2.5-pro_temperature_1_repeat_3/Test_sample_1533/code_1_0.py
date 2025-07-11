import math

def solve_geometry_problem():
    """
    This script provides a step-by-step solution to the geometry problem
    and prints the final result.
    """
    explanation = """
Here is the step-by-step derivation to find the ratio BM/MI in terms of the side lengths a, b, and c.

1.  **Identify Key Geometric Properties & The Trillium Theorem**
    - We are given a triangle ABC with side lengths a, b, c opposite to vertices A, B, C.
    - I is the incenter (intersection of angle bisectors).
    - M is the point where the angle bisector of ∠B (the line BI) intersects the circumcircle of triangle ABC.
    - A key result in geometry, known as the Incenter-Excenter Lemma or Trillium Theorem, states that M is equidistant from the vertices A, C, and the incenter I.
    - This gives us the crucial equality: MI = MA.

2.  **Simplify the Ratio**
    - The ratio we need to find is BM / MI.
    - Using the result from the Trillium Theorem (MI = MA), we can rewrite the ratio as:
      BM / MI = BM / MA

3.  **Apply the Law of Sines to Triangle ABM**
    - Let's analyze the triangle ABM.
    - The side opposite to angle AMB is AB = c.
    - The angles of triangle ABM are:
        - ∠ABM = B/2, because BM is the angle bisector of ∠ABC.
        - ∠AMB = ∠ACB = C, because both angles subtend the same arc AB on the circumcircle.
        - ∠MAB = 180° - (∠ABM + ∠AMB) = 180° - (B/2 + C).
          Since A + B + C = 180° in triangle ABC, we can rewrite this as:
          ∠MAB = (A + B + C) - (B/2 + C) = A + B/2.

    - Now, we apply the Law of Sines to triangle ABM:
      BM / sin(∠MAB) = MA / sin(∠ABM)
      BM / sin(A + B/2) = MA / sin(B/2)

4.  **Calculate the Ratio**
    - By rearranging the equation from the Law of Sines, we get our ratio:
      BM / MA = sin(A + B/2) / sin(B/2)

5.  **Convert the Trigonometric Expression to Side Lengths**
    - We have found that BM/MI = sin(A + B/2) / sin(B/2). The final step is to express this using side lengths a, b, and c.
    - Let's prove that sin(A + B/2) / sin(B/2) is equal to (a + c) / b.
    - First, let's express (a + c) / b using the Law of Sines on triangle ABC (where a/sin(A) = b/sin(B) = c/sin(C) = 2R):
      (a + c) / b = (2R*sin(A) + 2R*sin(C)) / (2R*sin(B)) = (sin(A) + sin(C)) / sin(B).
    - Using the sum-to-product and double-angle identities:
      (sin(A) + sin(C)) / sin(B) = [2*sin((A+C)/2)*cos((A-C)/2)] / [2*sin(B/2)*cos(B/2)].
    - Since A+B+C=180°, we have (A+C)/2 = 90° - B/2, so sin((A+C)/2) = cos(B/2).
      The expression simplifies to: [cos(B/2)*cos((A-C)/2)] / [sin(B/2)*cos(B/2)] = cos((A-C)/2) / sin(B/2).
    - Now let's simplify our ratio sin(A + B/2) / sin(B/2).
    - We use the identity A+B+C=180° to write A+B/2 = (A+B+C)/2 + (A-C)/2 = 90° + (A-C)/2.
    - Therefore, sin(A + B/2) = sin(90° + (A-C)/2) = cos((A-C)/2).
    - So, the ratio sin(A + B/2) / sin(B/2) also simplifies to cos((A-C)/2) / sin(B/2).
    - Since both BM/MI and (a+c)/b simplify to the same expression, they must be equal.

6.  **Final Result**
    - We have rigorously shown that the ratio is equivalent to the expression in terms of side lengths.
"""
    
    print(explanation)
    
    print("="*70)
    print("Final Answer")
    print("The ratio BM / MI can be expressed in terms of the side lengths a, b, and c as:")
    
    # The final equation is symbolic. We print its components as requested.
    numerator_part1 = "a"
    numerator_part2 = "c"
    denominator = "b"
    
    print(f"\n      {numerator_part1} + {numerator_part2}")
    print(f"Ratio = ———————")
    print(f"        {denominator}")

solve_geometry_problem()
<<< (a + c) / b >>>