import numpy as np

def solve():
    """
    Calculates the probability that the conditioned random walk never enters the set 
    of the four neighbors of the origin.
    """
    # Starting point
    x0 = (3000, 4000)

    # Calculate the Euclidean norm of the starting point
    norm_x0 = np.sqrt(x0[0]**2 + x0[1]**2)

    # Constants for the potential kernel formula
    pi = np.pi
    gamma = np.euler_gamma # Euler-Mascheroni constant
    log_8 = np.log(8)

    # Asymptotic formula for the potential kernel a(x) for large |x|
    # a(x) = (2/pi) * log|x| + (2*gamma + log(8))/pi
    term1 = (2 / pi) * np.log(norm_x0)
    term2 = (2 * gamma + log_8) / pi
    a_x0 = term1 + term2

    # The probability of the conditioned walk hitting the neighbors of the origin is ~ 1/a(x0)
    # The probability of never hitting is therefore 1 - 1/a(x0)
    prob = 1 - (1 / a_x0)

    print("Step 1: Calculate the norm of the starting point.")
    print(f"x_0 = {x0}")
    print(f"||x_0|| = {norm_x0:.0f}\n")

    print("Step 2: Calculate the value of the potential kernel a(x_0) using the asymptotic formula.")
    print(f"a(x_0) ≈ (2/π) * log(||x_0||) + (2γ + log(8))/π")
    print(f"a({x0}) ≈ (2/{pi:.4f}) * log({norm_x0:.0f}) + (2*{gamma:.4f} + {log_8:.4f})/{pi:.4f}")
    print(f"a({x0}) ≈ {term1:.4f} + {term2:.4f} = {a_x0:.4f}\n")

    print("Step 3: Calculate the final probability.")
    print(f"P(never hit N) ≈ 1 - 1/a(x_0)")
    print(f"P ≈ 1 - 1/{a_x0:.4f}")
    print(f"P ≈ {prob:.4f}\n")

    print(f"The approximate probability with two significant digits is {prob:.2f}.")

solve()