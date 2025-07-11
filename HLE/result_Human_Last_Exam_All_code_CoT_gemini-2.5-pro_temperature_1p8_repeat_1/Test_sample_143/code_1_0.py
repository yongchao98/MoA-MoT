import numpy as np

def check_convexity_counterexample():
    """
    This function demonstrates with a counter-example that for d=2, the set of
    achievable squared projection payoffs is not convex, even with non-negative data.
    """
    # Define the data points x_i for d=2, n=3
    x1 = np.array([1., 0.])
    x2 = np.array([0., 1.])
    x3 = np.array([1., 1.])

    # Define a function to compute the payoff vector for a given w
    def get_payoff(w):
        y1 = (x1 @ w)**2
        y2 = (x2 @ w)**2
        y3 = (x3 @ w)**2
        return np.array([y1, y2, y3])

    # 1. Take two points, A and B, from the payoff set Z
    w_A = np.array([1., 0.]) # Corresponds to theta = 0
    Y_A = get_payoff(w_A)

    w_B = np.array([0., 1.]) # Corresponds to theta = pi/2
    Y_B = get_payoff(w_B)

    # 2. Calculate their midpoint Y_M
    Y_M = 0.5 * Y_A + 0.5 * Y_B
    
    print("This code provides a counterexample for d=2.")
    print(f"Let x1 = {x1}, x2 = {x2}, x3 = {x3}.")
    print(f"For wA = {w_A}, the payoff vector YA = {Y_A}")
    print(f"For wB = {w_B}, the payoff vector YB = {Y_B}")
    print(f"The midpoint is YM = {Y_M}")
    print("\nTo check if YM is in the payoff set, we need to find w=(w1, w2) with w1^2+w2^2=1 such that:")
    print(f"(x1.w)^2 = w1^2 = {Y_M[0]}")
    print(f"(x2.w)^2 = w2^2 = {Y_M[1]}")
    print(f"(x3.w)^2 = (w1+w2)^2 = {Y_M[2]}")

    # 3. Prove that Y_M is not in the set Z
    # We are looking for w=(w1, w2) satisfying the system of equations.
    # w1^2 = 0.5 => w1 = +/- sqrt(0.5)
    # w2^2 = 0.5 => w2 = +/- sqrt(0.5)
    # (w1+w2)^2 = 1 => w1^2 + w2^2 + 2*w1*w2 = 1
    # Substituting the first two equations into the third:
    # 0.5 + 0.5 + 2*w1*w2 = 1 => 1 + 2*w1*w2 = 1 => w1*w2 = 0.
    # This contradicts w1 = +/- sqrt(0.5) and w2 = +/- sqrt(0.5).
    # Thus, no such w exists.

    print("\nSolving this system of equations leads to a contradiction:")
    print("From the first two equations, w1 and w2 must both be non-zero.")
    print("The third equation simplifies to w1*w2 = 0, which means at least one must be zero.")
    print("Since there is a contradiction, no such w exists.")
    print("This proves that the midpoint YM is not in the set of achievable payoffs.")
    print("Therefore, the set is not convex for d=2, even with non-negative data.")

check_convexity_counterexample()