import math

def solve_branching_walk():
    """
    Calculates the limit of the probability that site 0 is visited by infinitely
    many different particles in a branching random walk in random environment.
    """

    # Let's define the parameters of the problem.
    # The probability of a site being red is h.
    # The probability of a particle branching is also h.
    # We are interested in the limit as h -> 0.
    
    # We can set h to its limit value for the calculation.
    h = 0

    # Define jump probabilities for red and blue sites.
    p_L_red = 4/5
    p_R_red = 1/5
    p_L_blue = 1/5
    p_R_blue = 4/5

    # The behavior of the walk is determined by the average logarithmic drift.
    # Let rho = p_L / p_R be the ratio of left to right jump probabilities.
    # We calculate the expected value of log(rho).
    
    # At a red site, rho_red = (4/5) / (1/5) = 4.
    rho_red = p_L_red / p_R_red
    log_rho_red = math.log(rho_red)
    
    # At a blue site, rho_blue = (1/5) / (4/5) = 1/4.
    rho_blue = p_L_blue / p_R_blue
    log_rho_blue = math.log(rho_blue)

    # The expected value of log(rho) is E[log(rho)] = P(red)*log(rho_red) + P(blue)*log(rho_blue).
    # This equation is: E[log(rho)] = h * log(4) + (1-h) * log(1/4)
    # which simplifies to: E[log(rho)] = h * log(4) - (1-h) * log(4) = log(4) * (2*h - 1).
    
    # Now we take the limit of this expectation as h -> 0.
    expected_log_rho_limit = math.log(4) * (2 * h - 1)

    print("Step 1: Determine if the walk is transient or recurrent.")
    print("This depends on the sign of E[log(p_L/p_R)].")
    print(f"For a red site, p_L/p_R = ({p_L_red})/({p_R_red}) = {rho_red}. log(p_L/p_R) = {log_rho_red:.4f}")
    print(f"For a blue site, p_L/p_R = ({p_L_blue})/({p_R_blue}) = {rho_blue:.4f}. log(p_L/p_R) = {log_rho_blue:.4f}")
    print("\nStep 2: Calculate the expectation in the limit h -> 0.")
    print("E[log(rho)] is given by the expression: log(4) * (2*h - 1)")
    print(f"lim_{{h->0}} E[log(rho)] = log(4) * (2*0 - 1) = {expected_log_rho_limit:.4f}")
    
    print("\nStep 3: Interpret the result.")
    print(f"Since the expected log-ratio ({expected_log_rho_limit:.4f}) is negative, the random walk is transient to the right.")
    print("This means the entire particle cloud drifts towards +infinity.")
    print("As a result, any specific site, including site 0, will be visited only a finite number of times.")
    
    # The probability of infinitely many visits is therefore 0.
    final_result = 0

    print("\nFinal Conclusion:")
    print("The probability of site 0 being visited by infinitely many different particles is 0.")
    print("The final equation for the limit is:")
    print(f"lim_{{h->0}} P[site 0 is visited by infinitely many different particles] = {final_result}")

solve_branching_walk()