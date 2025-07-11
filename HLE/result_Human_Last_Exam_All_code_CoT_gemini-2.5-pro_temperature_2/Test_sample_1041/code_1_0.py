import math

def calculate_hitting_probability():
    """
    Calculates the probability for a 2D random walk to hit a target set
    before escaping a disk, using an approximation from potential theory.
    """
    
    # Define the parameters of the problem
    R = 1000.0  # Radius of the escape disk
    start_point = (0, 300)
    target_point_1 = (0, 0)
    target_point_2 = (2, 0)

    # Calculate the distance 'd' from the starting point to the target set.
    # Since the start point is far from the target set, we can use the distance
    # to one of the points as a close approximation.
    d = math.sqrt((start_point[0] - target_point_1[0])**2 + (start_point[1] - target_point_1[1])**2)
    
    # Calculate the effective radius 'epsilon_eff' of the target set.
    # For a target of two points separated by distance s, this is ~sqrt(s).
    s = math.sqrt((target_point_2[0] - target_point_1[0])**2 + (target_point_2[1] - target_point_1[1])**2)
    epsilon_eff = math.sqrt(s)
    
    # The probability `p` is given by the formula based on logarithmic potential:
    # p = (log(R) - log(d)) / (log(R) - log(epsilon_eff))
    
    # Calculate logarithms for each term
    log_R = math.log(R)
    log_d = math.log(d)
    log_epsilon_eff = math.log(epsilon_eff)
    
    # Calculate numerator and denominator
    numerator = log_R - log_d
    denominator = log_R - log_epsilon_eff
    
    # Calculate the final probability
    probability = numerator / denominator

    # Output the steps of the calculation as requested
    print("The hitting probability `p` is calculated using the formula:")
    print("p = (log(R) - log(d)) / (log(R) - log(epsilon_eff))\n")
    print("Where the values are:")
    print(f"R (escape radius)            = {R}")
    print(f"d (distance to target)       = {d:.3f}")
    print(f"epsilon_eff (effective radius) = {epsilon_eff:.3f}\n")
    
    print("Plugging these into the formula:")
    print(f"p = (log({R}) - log({d:.3f})) / (log({R}) - log({epsilon_eff:.3f}))")
    print(f"p = ({log_R:.3f} - {log_d:.3f}) / ({log_R:.3f} - {log_epsilon_eff:.3f})")
    print(f"p = {numerator:.3f} / {denominator:.3f}\n")

    print(f"The calculated probability is: {probability:.5f}")
    print(f"The probability rounded to three significant digits is: {probability:.3f}")


calculate_hitting_probability()