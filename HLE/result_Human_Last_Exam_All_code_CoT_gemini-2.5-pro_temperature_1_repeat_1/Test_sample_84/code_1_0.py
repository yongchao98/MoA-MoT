import math

def calculate_alpha():
    """
    This function calculates the value of alpha based on the problem's parameters.

    The logic is based on analyzing the required degree of a polynomial that
    is small on one interval and large on another. The degree 'd' is found
    to be proportional to 1/sqrt(epsilon), where epsilon is the scaled size
    of the gap between the intervals.
    """

    # The problem defines two sets of points based on integer values up to n^2 and n^10.
    # The length of the first interval of constraints is on the order of n^2.
    # The length of the second interval of constraints is on the order of n^10.
    n_exponent_S1 = 2
    n_exponent_S2 = 10

    # We map the larger interval (related to n^10) to [-1, 1].
    # The smaller interval (related to n^2) is mapped to a region just outside [-1, 1].
    # The distance from the end of the mapped smaller interval to the nearest point in [-1, 1]
    # is denoted by epsilon. This distance scales inversely with the length of the larger interval.
    # epsilon is proportional to 1 / n^n_exponent_S2
    epsilon_exponent = n_exponent_S2

    # The degree 'd' of the Chebyshev polynomial required to achieve the necessary
    # amplification over the small distance epsilon is inversely proportional to sqrt(epsilon).
    # d ~ 1 / sqrt(epsilon)
    # d ~ 1 / sqrt(1 / n^epsilon_exponent)
    # d ~ 1 / (n^(-epsilon_exponent/2))
    # d ~ n^(epsilon_exponent/2)
    # So, the exponent alpha is half of the exponent of the larger interval's length.
    
    alpha = n_exponent_S2 / 2

    print(f"The first set of constraints is on points up to n^{n_exponent_S1}.")
    print(f"The second set of constraints is on points up to n^{n_exponent_S2}.")
    print(f"The asymptotic growth rate of the polynomial degree d_n is of the form Theta(n^alpha).")
    print(f"The value of alpha is determined by the exponent of the larger interval.")
    print(f"alpha = n_exponent_S2 / 2")
    print(f"alpha = {n_exponent_S2} / 2")
    print(f"alpha = {int(alpha)}")
    
    return int(alpha)

# Run the calculation and print the final answer.
final_alpha = calculate_alpha()