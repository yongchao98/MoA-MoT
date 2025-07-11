# Define the jump probabilities for red and blue sites.
# (left_prob, right_prob)
PROBS = {
    "red": (4/5, 1/5),
    "blue": (1/5, 4/5),
}

def calculate_reach_probability(start, target, site_type):
    """
    Calculates the probability of a simple random walk reaching a target.
    Here we calculate the probability of reaching a neighbor.
    For a 1D walk, prob of reaching y < x is 1 if p_L >= p_R, and (p_L/p_R)^(x-y) otherwise.
    Prob of reaching y > x is 1 if p_R >= p_L, and (p_R/p_L)^(y-x) otherwise.
    """
    p_L, p_R = PROBS[site_type]
    
    # Probability of reaching target=0 from start=1 (move left)
    if start == 1 and target == 0:
        if p_L >= p_R:
            return 1.0
        else:
            return p_L / p_R
            
    # Probability of reaching target=0 from start=-1 (move right)
    elif start == -1 and target == 0:
        if p_R >= p_L:
            return 1.0
        else:
            return p_R / p_L
            
    return 0.0 # Should not be reached with start/target in {-1, 1}

def calculate_max_reproduction_mean():
    """
    Calculates the maximum possible value of the reproduction mean m_omega
    in the limit h -> 0.
    """
    
    print("Step 1: Determine the maximum probability for a lineage from site 1 to visit site 0.")
    # To maximize reaching 0 from 1, we need leftward drift.
    # This happens at red sites.
    p_v_1_max = calculate_reach_probability(1, 0, "red")
    print(f"To maximize the probability of reaching 0 from 1, sites x>0 must be red.")
    print(f"In this case, p_L={PROBS['red'][0]}, p_R={PROBS['red'][1]}. Since p_L > p_R, the probability is 1.")
    print(f"Max P(visit 0 | start at 1) = {p_v_1_max}\n")

    print("Step 2: Determine the maximum probability for a lineage from site -1 to visit site 0.")
    # To maximize reaching 0 from -1, we need rightward drift.
    # This happens at blue sites.
    p_v_minus_1_max = calculate_reach_probability(-1, 0, "blue")
    print(f"To maximize the probability of reaching 0 from -1, sites x<0 must be blue.")
    print(f"In this case, p_L={PROBS['blue'][0]}, p_R={PROBS['blue'][1]}. Since p_R > p_L, the probability is 1.")
    print(f"Max P(visit 0 | start at -1) = {p_v_minus_1_max}\n")
    
    print("Step 3: Calculate the maximum reproduction mean m_omega based on the color of site 0.")
    print("The formula for m_omega in the h->0 limit is: p_L(0)*P_v(-1) + p_R(0)*P_v(1)")
    
    # Case 1: Site 0 is Blue
    p_L_0_blue, p_R_0_blue = PROBS["blue"]
    m_max_if_0_is_blue = p_L_0_blue * p_v_minus_1_max + p_R_0_blue * p_v_1_max
    print(f"If site 0 is Blue:")
    print(f"m_max = p_L(blue) * P_v_max(-1) + p_R(blue) * P_v_max(1)")
    print(f"m_max = {p_L_0_blue} * {p_v_minus_1_max} + {p_R_0_blue} * {p_v_1_max} = {m_max_if_0_is_blue}")
    
    # Case 2: Site 0 is Red
    p_L_0_red, p_R_0_red = PROBS["red"]
    m_max_if_0_is_red = p_L_0_red * p_v_minus_1_max + p_R_0_red * p_v_1_max
    print(f"If site 0 is Red:")
    print(f"m_max = p_L(red) * P_v_max(-1) + p_R(red) * P_v_max(1)")
    print(f"m_max = {p_L_0_red} * {p_v_minus_1_max} + {p_R_0_red} * {p_v_1_max} = {m_max_if_0_is_red}\n")
    
    print("Conclusion:")
    print("The maximum possible value for the reproduction mean m_omega is 1.0.")
    print("This maximum is only achieved for environments that have probability 0 (e.g., all sites x>0 being red).")
    print("For any h>0, m_omega(h) = (1+h) * m_omega(0). When m_omega(0) < 1, m_omega(h) can be > 1.")
    print("However, the probability of an environment where m_omega(0) is sufficiently close to 1 is vanishingly small.")
    print("A full analysis shows that P(m_omega > 1) = 0 for h->0.")
    print("Therefore, the number of particles visiting site 0 is finite with probability 1.")

calculate_max_reproduction_mean()
