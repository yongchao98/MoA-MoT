import math
from fractions import Fraction

def solve_betting_problem():
    """
    Calculates and explains the solution to the bike racing betting problem.
    """
    # Define probabilities and odds using Fraction for precision
    p_true = [Fraction(1, 2), Fraction(1, 4), Fraction(1, 8), Fraction(1, 8)]
    q_incorrect = [Fraction(1, 4), Fraction(1, 2), Fraction(1, 8), Fraction(1, 8)]
    odds = [4, 3, 7, 7]

    # Use string representations for clear printing
    p_str = [str(f) for f in p_true]
    q_str = [str(f) for f in q_incorrect]

    print("--- Calculating the Optimal Growth Rate W* ---")
    print("The optimal strategy is b_i = p_i. The growth rate W* = Σ p_i * ln(p_i * o_i).")
    print(f"W* = ({p_str[0]})*ln(({p_str[0]})*4) + ({p_str[1]})*ln(({p_str[1]})*3) + ({p_str[2]})*ln(({p_str[2]})*7) + ({p_str[3]})*ln(({p_str[3]})*7)")
    print(f"W* = ({p_str[0]})*ln({p_true[0]*odds[0]}) + ({p_str[1]})*ln({p_true[1]*odds[1]}) + ({p_str[2]})*ln({p_true[2]*odds[2]}) + ({p_str[3]})*ln({p_true[3]*odds[3]})")
    print("W* = (1/2)*ln(2) + (1/4)*ln(3/4) + (1/8)*ln(7/8) + (1/8)*ln(7/8)")
    print("Simplified W* = (1/2)*ln(2) + (1/4)*ln(3/4) + (1/4)*ln(7/8)")
    print("\n" + "="*50 + "\n")

    print("--- Calculating the Achieved Growth Rate W ---")
    print("Your strategy is b_i = q_i. The achieved growth rate W = Σ p_i * ln(q_i * o_i).")
    print(f"W = ({p_str[0]})*ln(({q_str[0]})*4) + ({p_str[1]})*ln(({q_str[1]})*3) + ({p_str[2]})*ln(({q_str[2]})*7) + ({p_str[3]})*ln(({q_str[3]})*7)")
    print(f"W = ({p_str[0]})*ln({q_incorrect[0]*odds[0]}) + ({p_str[1]})*ln({q_incorrect[1]*odds[1]}) + ({p_str[2]})*ln({q_incorrect[2]*odds[2]}) + ({p_str[3]})*ln({q_incorrect[3]*odds[3]})")
    print("W = (1/2)*ln(1) + (1/4)*ln(3/2) + (1/8)*ln(7/8) + (1/8)*ln(7/8)")
    print("Since ln(1) = 0, we get:")
    print("W = (1/4)*ln(3/2) + (1/4)*ln(7/8)")
    print("\n" + "="*50 + "\n")

    print("--- Calculating the Decrease in Growth Rate ΔW ---")
    print("The decrease ΔW = W* - W is the KL-divergence D(p || q) = Σ p_i * ln(p_i / q_i).")
    p_over_q = [p_true[i] / q_incorrect[i] for i in range(4)]
    p_over_q_str = [str(f) for f in p_over_q]
    print(f"ΔW = ({p_str[0]})*ln({p_over_q_str[0]}) + ({p_str[1]})*ln({p_over_q_str[1]}) + ({p_str[2]})*ln({p_over_q_str[2]}) + ({p_str[3]})*ln({p_over_q_str[3]})")
    print("ΔW = (1/2)*ln(2) + (1/4)*ln(1/2) + (1/8)*ln(1) + (1/8)*ln(1)")
    print("Since ln(1/2) = -ln(2) and ln(1) = 0:")
    print("ΔW = (1/2)*ln(2) - (1/4)*ln(2)")
    print("ΔW = (1/4)*ln(2)")

if __name__ == '__main__':
    solve_betting_problem()