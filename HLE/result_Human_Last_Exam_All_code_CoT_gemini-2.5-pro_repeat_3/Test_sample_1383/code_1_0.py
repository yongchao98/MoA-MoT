import math
from fractions import Fraction

def solve_betting_problem():
    """
    Calculates and explains the achieved growth rate (W) and the decrease in
    growth rate (Delta W) for the given betting problem.
    """
    # Define the true probabilities, believed probabilities, and odds
    p = [Fraction(1, 2), Fraction(1, 4), Fraction(1, 8), Fraction(1, 8)]
    q = [Fraction(1, 4), Fraction(1, 2), Fraction(1, 8), Fraction(1, 8)]
    # Odds are O-for-1, so the total return is O times the bet.
    O = [4, 3, 7, 7]

    print("--- Part 1: Calculating the Achieved Growth Rate W ---")
    print("You bet according to your incorrect probabilities q, so your betting fractions are f_i = q_i.")
    print("The achieved growth rate W is calculated with the true probabilities p_i.")
    print("Formula: W = p1*log(q1*O1) + p2*log(q2*O2) + p3*log(q3*O3) + p4*log(q4*O4)\n")

    print("Substituting the given values:")
    w_equation = (
        f"W = ({p[0]})*log(({q[0]})*{O[0]}) + ({p[1]})*log(({q[1]})*{O[1]}) + "
        f"({p[2]})*log(({q[2]})*{O[2]}) + ({p[3]})*log(({q[3]})*{O[3]})"
    )
    print(w_equation)

    print("\nSimplifying the terms inside the logarithms:")
    w_simplified_terms = [f"({p[i]})*log({q[i]*O[i]})" for i in range(4)]
    print("W = " + " + ".join(w_simplified_terms))

    print("\nAfter further simplification (since log(1)=0 and combining the last two terms):")
    final_w_expr = "W = (1/4) * log(3/2) + (1/4) * log(7/8)"
    print(final_w_expr)

    print("\n" + "="*70)

    print("--- Part 2: Calculating the Decrease in Growth Rate (Delta W) ---")
    print("The decrease, Delta W, is the difference between the optimal rate W* and the achieved rate W.")
    print("This is equivalent to the Kullback-Leibler divergence D_KL(p || q).")
    print("Formula: Delta W = p1*log(p1/q1) + p2*log(p2/q2) + p3*log(p3/q3) + p4*log(p4/q4)\n")
    
    print("Substituting the given values:")
    delta_w_equation_parts = []
    for i in range(4):
        delta_w_equation_parts.append(f"({p[i]})*log(({p[i]})/({q[i]}))")
    print("Delta W = " + " + ".join(delta_w_equation_parts))

    print("\nSimplifying the terms inside the logarithms:")
    delta_w_simplified_terms = [f"({p[i]})*log({p[i]/q[i]})" for i in range(4)]
    print("Delta W = " + " + ".join(delta_w_simplified_terms))

    print("\nAfter further simplification (since log(1)=0 and log(1/2)=-log(2)):")
    final_delta_w_expr = "Delta W = (1/4) * log(2)"
    print(final_delta_w_expr)

# Execute the function to print the solution
solve_betting_problem()