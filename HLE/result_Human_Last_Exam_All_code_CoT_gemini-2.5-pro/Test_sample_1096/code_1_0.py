import math

def solve_optimal_steps():
    """
    Calculates the optimal step sizes (gamma_1, gamma_2) for a two-step
    gradient descent method on an M-smooth and mu-strongly convex function.
    """
    # The problem specifies using m=1 (mu=1) and M=kappa (the condition number).
    # We will use M=10.0 as an example value for the condition number.
    M = 10.0
    mu = 1.0

    print(f"For an M-smooth and mu-strongly convex function with M = {M} and mu = {mu}:")

    # The optimal parameters are derived from the roots of the scaled Chebyshev
    # polynomial of degree 2. The step sizes gamma_1 and gamma_2 are the
    # reciprocals of these roots.
    # The formula for the roots is: r = (M+mu)/2 +- (M-mu)/(2*sqrt(2))

    sqrt2 = math.sqrt(2)
    term_sum = (M + mu) / 2
    term_diff = (M - mu) / (2 * sqrt2)

    # Calculate the two roots
    root1 = term_sum - term_diff
    root2 = term_sum + term_diff

    # The optimal step sizes are the reciprocals of these roots.
    gamma1 = 1 / root1
    gamma2 = 1 / root2

    # The final answer is the pair (gamma1, gamma2). The order does not matter.
    # As requested, we show the equation with all the numbers.
    print("\nThe best choice for the pair (gamma_1, gamma_2) is calculated as follows:")
    print("---------------------------------------------------------------------")

    print("\nCalculation for the first step size (let's call it gamma_a):")
    print(f"gamma_a = 1 / ( ({M} + {mu})/2 - ({M} - {mu})/(2 * sqrt(2)) )")
    print(f"gamma_a = 1 / ( {term_sum:.4f} - {term_diff:.4f} )")
    print(f"gamma_a = 1 / {root1:.4f}")
    print(f"gamma_a = {gamma1:.4f}")

    print("\nCalculation for the second step size (let's call it gamma_b):")
    print(f"gamma_b = 1 / ( ({M} + {mu})/2 + ({M} - {mu})/(2 * sqrt(2)) )")
    print(f"gamma_b = 1 / ( {term_sum:.4f} + {term_diff:.4f} )")
    print(f"gamma_b = 1 / {root2:.4f}")
    print(f"gamma_b = {gamma2:.4f}")

    print("\n---------------------------------------------------------------------")
    print(f"The optimal pair of step sizes (gamma_1, gamma_2) is ({gamma1:.4f}, {gamma2:.4f}).")


if __name__ == '__main__':
    solve_optimal_steps()
