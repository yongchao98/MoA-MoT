import math
from fractions import Fraction

def solve_betting_problem():
    """
    Calculates and displays the results for the bike racing betting problem.
    """
    # True probabilities
    p = [Fraction(1, 2), Fraction(1, 4), Fraction(1, 8), Fraction(1, 8)]
    # Incorrectly believed probabilities
    q = [Fraction(1, 4), Fraction(1, 2), Fraction(1, 8), Fraction(1, 8)]
    # Decimal odds (X-for-1)
    d = [4, 3, 7, 7]

    print("This script calculates the expected logarithmic growth rates for a horse race betting scenario.")
    print("-" * 80)

    # --- 1. Calculate the achieved doubling rate W ---
    # W = Σ pᵢ * log(qᵢ * dᵢ)
    # Term 1: p₁*log(q₁*d₁) = (1/2)*log(1/4 * 4) = (1/2)*log(1) = 0
    # Term 2: p₂*log(q₂*d₂) = (1/4)*log(1/2 * 3) = (1/4)*log(3/2)
    # Term 3: p₃*log(q₃*d₃) = (1/8)*log(1/8 * 7) = (1/8)*log(7/8)
    # Term 4: p₄*log(q₄*d₄) = (1/8)*log(1/8 * 7) = (1/8)*log(7/8)
    # Combining terms 3 and 4: (1/8 + 1/8)*log(7/8) = (1/4)*log(7/8)
    # Final W = (1/4)*log(3/2) + (1/4)*log(7/8)
    
    w_coeff_1 = p[1]
    w_arg_1 = q[1] * d[1]
    w_coeff_2 = p[2] + p[3]
    w_arg_2 = q[2] * d[2]
    
    print("The achieved doubling rate W (using incorrect probabilities) is:")
    print(f"W = ({p[0]})*log(({q[0]})*({d[0]})) + ({p[1]})*log(({q[1]})*({d[1]})) + ({p[2]})*log(({q[2]})*({d[2]})) + ({p[3]})*log(({q[3]})*({d[3]}))")
    print("Simplified, this is:")
    print(f"W = ({w_coeff_1}) * log({w_arg_1}) + ({w_coeff_2}) * log({w_arg_2})")
    print("-" * 80)
    
    # --- 2. Calculate the optimal doubling rate W* ---
    # W* = Σ pᵢ * log(pᵢ * dᵢ)
    # Term 1: p₁*log(p₁*d₁) = (1/2)*log(1/2 * 4) = (1/2)*log(2)
    # Term 2: p₂*log(p₂*d₂) = (1/4)*log(1/4 * 3) = (1/4)*log(3/4)
    # Term 3 & 4: (1/4)*log(7/8)
    
    w_star_coeff_1 = p[0]
    w_star_arg_1 = p[0] * d[0]
    w_star_coeff_2 = p[1]
    w_star_arg_2 = p[1] * d[1]
    w_star_coeff_3 = p[2] + p[3]
    w_star_arg_3 = p[2] * d[2]
    
    print("The optimal doubling rate W* (using true probabilities) is:")
    print(f"W* = ({p[0]})*log(({p[0]})*({d[0]})) + ({p[1]})*log(({p[1]})*({d[1]})) + ({p[2]})*log(({p[2]})*({d[2]})) + ({p[3]})*log(({p[3]})*({d[3]}))")
    print("Simplified, this is:")
    print(f"W* = ({w_star_coeff_1}) * log({w_star_arg_1}) + ({w_star_coeff_2}) * log({w_star_arg_2}) + ({w_star_coeff_3}) * log({w_star_arg_3})")
    print("-" * 80)

    # --- 3. Calculate the decrease in rate ΔW ---
    # ΔW = W* - W = D_KL(p || q) = Σ pᵢ * log(pᵢ/qᵢ)
    # Term 1: p₁*log(p₁/q₁) = (1/2)*log((1/2)/(1/4)) = (1/2)*log(2)
    # Term 2: p₂*log(p₂/q₂) = (1/4)*log((1/4)/(1/2)) = (1/4)*log(1/2) = -(1/4)*log(2)
    # Term 3: p₃*log(p₃/q₃) = (1/8)*log((1/8)/(1/8)) = (1/8)*log(1) = 0
    # Term 4: p₄*log(p₄/q₄) = (1/8)*log((1/8)/(1/8)) = (1/8)*log(1) = 0
    # Final ΔW = (1/2)*log(2) - (1/4)*log(2) = (1/4)*log(2)
    
    delta_w_coeff = Fraction(1, 4)
    delta_w_arg = 2
    
    print("The decrease in the doubling rate ΔW = W* - W is:")
    print(f"ΔW = ({delta_w_coeff}) * log({delta_w_arg})")
    print("-" * 80)

solve_betting_problem()