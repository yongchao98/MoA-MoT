import math

def solve_random_walk_probability():
    """
    Calculates the approximate probability for the described random walk problem.
    """
    # Define constants
    pi = math.pi
    # The Euler-Mascheroni constant
    gamma = 0.5772156649015328

    # Define norms of the start and target points
    x0_norm = 1
    y_norm = 3600

    # The problem is to find the probability that a 2D simple random walk,
    # starting at x0=(0,1) and conditioned to never enter the origin,
    # eventually hits the set of four neighbours of y=(3600,0).

    # This probability P can be approximated using the potential kernel a(x).
    # The approximation is P ≈ 4 * a(x0) / a(y).

    # Step 1: Calculate the potential kernel at the starting point x0=(0,1).
    # For a neighbor of the origin, the exact value is known.
    # Formula: a(x0) = 4/π - 1
    a_x0 = 4 / pi - 1

    # Step 2: Calculate the potential kernel at the target point y=(3600,0).
    # For large |y|, we use the asymptotic formula:
    # a(y) ≈ (2/π)ln|y| + (2γ + ln(8))/π
    ln_y = math.log(y_norm)
    ln_8 = math.log(8)
    # The constant term C in the asymptotic expansion
    C = (2 * gamma + ln_8) / pi
    a_y = (2 / pi) * ln_y + C

    # Step 3: Calculate the final probability using the approximation.
    prob = 4 * a_x0 / a_y

    # Output the explanation and the step-by-step calculation.
    print("The probability P is approximated by the formula: P ≈ 4 * a(x0) / a(y)")
    print("-" * 50)
    print("Step 1: Calculation of a(x0) for x0=(0,1)")
    print(f"a(x0) = 4 / π - 1")
    print(f"a(x0) = 4 / {pi:.6f} - 1 = {a_x0:.6f}")
    print("-" * 50)
    print("Step 2: Calculation of a(y) for y=(3600,0)")
    print(f"a(y) ≈ (2/π) * ln(|y|) + (2γ + ln(8))/π")
    print(f"a(y) ≈ (2/{pi:.6f}) * ln({y_norm}) + (2*{gamma:.6f} + {ln_8:.6f})/{pi:.6f}")
    print(f"a(y) ≈ {a_y:.6f}")
    print("-" * 50)
    print("Step 3: Final probability calculation")
    print(f"P ≈ 4 * a(x0) / a(y) = 4 * {a_x0:.6f} / {a_y:.6f}")
    print(f"P ≈ {prob:.6f}")
    print("-" * 50)
    # Print the final answer with two significant digits
    print(f"The approximate probability is {prob:.2g}.")

solve_random_walk_probability()