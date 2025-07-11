from fractions import Fraction

def solve_stingray_problem():
    """
    Solves the stingray duo problem and prints the results.
    """
    d = 5
    q = 4
    e1 = 3
    e2 = 2

    # Part (a): Is the pair irreducible?
    # A (3,2)-stingray duo is not always irreducible. It can be reducible if, for example,
    # the fixed space of g1 intersects the fixed space of g2 non-trivially.
    answer_a = "No"

    # Part (b): Causes of reducibility.
    # The standard theorem states that for a stingray duo with e1+e2=d, reducibility
    # occurs if and only if one of the three listed conditions holds.
    answer_b = "{(1), (2), (3)}"

    # Part (c): Calculate the proportion of irreducible (3,2)-stingray duos.
    # This is interpreted as the probability that a random (3,2)-stingray duo is irreducible.
    # P(irr) = 1 - P(reducible)
    # The events for reducibility are mutually exclusive.
    # P(reducible) = P(F1 intersect F2 != 0) + P(U1 = F2) + P(U2 = F1)

    print("(a) {}".format(answer_a))
    print("(b) {}".format(answer_b))
    print("(c) The proportion of irreducible (3,2)-stingray duos among all such duos is calculated as follows:")
    print("Let P_irr be the proportion. P_irr = 1 - P_red, where P_red is the proportion of reducible duos.")
    print("The three causes for reducibility are mutually exclusive events.")
    print("P_red = P(F1 \u2229 F2 \u2260 {0}) + P(U1 = F2) + P(U2 = F1)")
    print("-" * 20)

    # P(F1 intersect F2 != 0) = 1/q
    p_f_intersect = Fraction(1, q)
    print(f"P(F1 \u2229 F2 \u2260 {{0}}) = 1/q = 1/{q}")

    # P(U2 = F1) = 1 / q^(e1 * (d-e1))
    exp1 = e1 * (d - e1)
    p_u2_f1_den = q**exp1
    p_u2_f1 = Fraction(1, p_u2_f1_den)
    print(f"P(U2 = F1) = 1 / q^(e1 * (d-e1)) = 1 / {q}^({e1}*({d}-{e1})) = 1 / {q}^{exp1} = 1/{p_u2_f1_den}")

    # P(U1 = F2) = 1 / q^(e2 * (d-e2))
    exp2 = e2 * (d - e2)
    p_u1_f2_den = q**exp2
    p_u1_f2 = Fraction(1, p_u1_f2_den)
    print(f"P(U1 = F2) = 1 / q^(e2 * (d-e2)) = 1 / {q}^({e2}*({d}-{e2})) = 1 / {q}^{exp2} = 1/{p_u1_f2_den}")
    print("-" * 20)

    # P_red = p_f_intersect + p_u2_f1 + p_u1_f2
    p_red = p_f_intersect + p_u2_f1 + p_u1_f2
    print(f"P_red = {p_f_intersect.numerator}/{p_f_intersect.denominator} + {p_u2_f1.numerator}/{p_u2_f1.denominator} + {p_u1_f2.numerator}/{p_u1_f2.denominator}")
    common_den = p_red.denominator
    num1 = p_f_intersect.numerator * (common_den // p_f_intersect.denominator)
    num2 = p_u2_f1.numerator * (common_den // p_u2_f1.denominator)
    num3 = p_u1_f2.numerator * (common_den // p_u1_f2.denominator)
    print(f"P_red = {num1}/{common_den} + {num2}/{common_den} + {num3}/{common_den} = {p_red.numerator}/{p_red.denominator}")
    
    # P_irr = 1 - P_red
    p_irr = 1 - p_red
    print(f"P_irr = 1 - {p_red.numerator}/{p_red.denominator} = {p_irr.numerator}/{p_irr.denominator}")
    
    final_answer_c = f"{p_irr.numerator}/{p_irr.denominator}"
    print("\nFinal Answer for (c):")
    print(final_answer_c)

solve_stingray_problem()