import collections

def solve_price():
    """
    Calculates the value of p based on the given probability distribution of effort and output.
    """
    # The given joint probability distribution P(e, y)
    # Format: { (effort, output): probability }
    prob_dist = {
        (22, 132): 0.4375,
        (22, 44): 0.0625,
        (10, 60): 0.0625,
        (10, 20): 0.4375,
    }

    # Calculate the expected output E[y]
    # E[y] = sum(y * P(e, y)) for all outcomes
    expected_y = 0.0
    for (e, y), prob in prob_dist.items():
        expected_y += y * prob
    
    # Calculate the expected squared effort E[e^2]
    # E[e^2] = sum(e^2 * P(e, y)) for all outcomes
    expected_e_squared = 0.0
    for (e, y), prob in prob_dist.items():
        expected_e_squared += (e**2) * prob

    # The price p is determined by the firm's profit maximization, which leads to
    # the condition p = E[e^2] / E[y]
    if expected_y == 0:
        p = float('inf') # Avoid division by zero, though not expected in this problem
    else:
        p = expected_e_squared / expected_y
        
    print(f"The calculated expected squared effort E[e^2] is: {expected_e_squared}")
    print(f"The calculated expected output E[y] is: {expected_y}")
    print("\nThe value of p is derived from the firm's optimality condition: p = E[e^2] / E[y]")
    print(f"p = {expected_e_squared} / {expected_y} = {p}")

solve_price()