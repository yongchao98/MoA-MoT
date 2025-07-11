import math
from fractions import Fraction

def solve():
    """
    Solves the betting problem based on the Kelly criterion for a superfair race.
    """
    # True probabilities
    p = {
        1: Fraction(1, 2),
        2: Fraction(1, 4),
        3: Fraction(1, 8),
        4: Fraction(1, 8)
    }

    # Incorrectly believed probabilities
    q = {
        1: Fraction(1, 4),
        2: Fraction(1, 2),
        3: Fraction(1, 8),
        4: Fraction(1, 8)
    }

    # Payout odds (e.g., 4-for-1 means o_1 = 4)
    o = {
        1: 4,
        2: 3,
        3: 7,
        4: 7
    }

    # --- Part 1: Optimal expected logarithmic growth rate W* ---
    # Optimal strategy is to bet fraction f_i = p_i
    # The return if bike i wins is p_i * o_i
    
    # Calculate returns for optimal strategy
    po = {i: p[i] * o[i] for i in p}
    
    # Expression for W*
    w_star_terms = []
    for i in sorted(p.keys()):
        w_star_terms.append(f"({p[i]})*log(({p[i]})*({o[i]}))")
    
    # Simplify terms for display
    w_star_display_terms = []
    for i in sorted(p.keys()):
        w_star_display_terms.append(f"{p[i]}*log({po[i]})")

    w_star_str = " + ".join(w_star_display_terms)
    
    # --- Part 2: Achieved doubling rate W ---
    # Strategy based on incorrect probabilities: f'_i = q_i
    # The return if bike i wins is q_i * o_i
    
    # Calculate returns for incorrect strategy
    qo = {i: q[i] * o[i] for i in q}

    # Expression for W, expectation taken over true probabilities p
    w_terms = []
    for i in sorted(p.keys()):
        w_terms.append(f"({p[i]})*log(({q[i]})*({o[i]}))")
        
    # Simplify terms for display
    w_display_terms = []
    for i in sorted(p.keys()):
        w_display_terms.append(f"{p[i]}*log({qo[i]})")

    w_str = " + ".join(w_display_terms)

    print("The doubling rate W you will achieve is:")
    print(f"W = {w_str}")
    
    # --- Part 3: Decrease in doubling rate delta_W = W* - W ---
    # delta_W = sum(p_i * log(p_i*o_i)) - sum(p_i * log(q_i*o_i))
    # delta_W = sum(p_i * log((p_i*o_i) / (q_i*o_i)))
    # delta_W = sum(p_i * log(p_i/q_i)) which is the KL divergence D_KL(p||q)

    kl_divergence_terms = []
    total_val = 0
    for i in sorted(p.keys()):
        # Avoid log(0) in case a probability is zero
        if p[i] > 0 and q[i] > 0:
            kl_divergence_terms.append(f"{p[i]}*log({p[i]}/{q[i]})")
            total_val += p[i] * math.log(p[i]/q[i])
        
    # Manually simplify the expression for a cleaner look
    # delta_W = (1/2)*log(2) + (1/4)*log(1/2) = (1/2)*log(2) - (1/4)*log(2) = (1/4)*log(2)
    delta_w_str = f"{p[1]}*log({p[1]/q[1]}) + {p[2]}*log({p[2]/q[2]})" # Only p1 and p2 are different
    delta_w_simplified_str = "1/4*log(2)"

    print("\nThe decrease in your doubling rate is:")
    print(f"ΔW = W* - W = {delta_w_simplified_str}")
    
solve()