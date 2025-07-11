import math

def solve():
    """
    Calculates the limit of the probability that site 0 is visited
    by infinitely many particles as h -> 0.
    """
    # The problem is to calculate lim_{h->0} P(site 0 is visited by infinitely many particles).
    # This probability is 0 if the expected total number of visits to site 0 is finite.
    # The expected total number of visits converges if (1+h) * 2 * sqrt(p*q) < 1.
    # We evaluate this expression in the limit h=0.

    h = 0

    # Probability of a site being red
    p_red = h
    # Probability of a site being blue
    p_blue = 1 - h

    # Jump probabilities at a red site
    p_left_red = 4/5
    p_right_red = 1/5

    # Jump probabilities at a blue site
    p_left_blue = 1/5
    p_right_blue = 4/5

    # Averaged jump probabilities
    q = p_red * p_left_red + p_blue * p_left_blue
    p = p_red * p_right_red + p_blue * p_right_blue

    # Mean number of offspring factor (each particle becomes 1+h on average)
    growth_factor = 1 + h

    # Convergence criterion
    convergence_param = growth_factor * 2 * math.sqrt(p * q)

    print("Step-by-step calculation for h=0:")
    print(f"h = {h}")
    print(f"Averaged probability of jumping right, p = (1-h)*4/5 + h*1/5 = {p}")
    print(f"Averaged probability of jumping left, q = (1-h)*1/5 + h*4/5 = {q}")
    print(f"Population growth factor per step, 1+h = {growth_factor}")
    print("\nFinal equation for the convergence criterion:")
    print(f"({growth_factor}) * 2 * sqrt(({p}) * ({q})) = {convergence_param}")

    # Since the convergence parameter is less than 1, the expected number of visits is finite.
    # This implies the probability of infinitely many visits is 0.
    limit_probability = 0
    print(f"\nSince {convergence_param:.2f} < 1, the expected number of total visits to site 0 is finite.")
    print("This implies the probability of infinitely many distinct particles visiting site 0 is 0.")
    print(f"The limit is: {limit_probability}")


solve()