import math

def solve():
    """
    Solves the committee election problem by calculating the satisfaction for W_1 and W_2
    and then finding their ratio.
    """
    
    # Define problem parameters based on the description
    num_voters = 10
    committee_size = 20
    group_N_size = 8
    
    # --- Step 1: Calculate s(N, W_1) for the core committee ---
    
    # The principle of proportionality suggests that a "core" committee allocates seats
    # based on group size.
    # Proportional share for group N = k * |N| / n
    k_N_core = committee_size * group_N_size / num_voters
    
    # The satisfaction for group N is s(N,W) = 7*w_c + k_N, where w_c is the number of
    # commonly approved candidates {1..8} selected for N's seats.
    # To find W_1, which gives the lowest satisfaction, we must minimize s(N,W).
    # This is achieved by minimizing w_c. The minimum possible value is w_c = 0.
    min_w_c_core = 0
    s_N_W1 = 7 * min_w_c_core + k_N_core
    
    # --- Step 2: Calculate s(N, W_2) for the EJR committee ---
    
    # To find W_2, we must find an EJR committee that minimizes s(N,W).
    # The satisfaction formula is s(N,W) = 7*w_c + k_N.
    # We need to find the minimal feasible values for w_c and k_N under EJR constraints.
    
    # The EJR requirement for group N (|A(i) intersect W| >= 8 for some i in N)
    # implies that w_c must be at least 6.
    min_w_c_ejr = 6
    
    # The EJR requirements for groups {9} and {10} set lower bounds on their seat counts.
    # To minimize k_N (number of seats for group N), we maximize the seats for {9} and {10}.
    # They each approve 4 candidates, so their max seats are 4 each.
    # min k_N = total_seats - max_k_9 - max_k_10
    min_k_N_ejr = committee_size - 4 - 4
    
    # This minimal configuration (w_c=6, k_N=12) is feasible under EJR.
    s_N_W2 = 7 * min_w_c_ejr + min_k_N_ejr
    
    # --- Step 3: Calculate and print the ratio ---
    
    # Ensure results are integers for the final output
    numerator = int(s_N_W1)
    denominator = int(s_N_W2)
    
    # Simplify the fraction by dividing by the greatest common divisor
    common_divisor = math.gcd(numerator, denominator)
    
    print(f"The lowest satisfaction for group N with a core committee is s(N,W_1) = {numerator}.")
    print(f"The lowest satisfaction for group N with an EJR committee is s(N,W_2) = {denominator}.")
    print(f"The final equation is: {numerator} / {denominator} = {numerator//common_divisor} / {denominator//common_divisor}")

solve()