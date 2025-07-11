import math

def solve_optimal_steps():
    """
    Calculates the optimal step sizes (gamma_1, gamma_2) for a two-step
    gradient descent method.
    """
    # The problem assumes mu=1 and M=kappa (the condition number).
    # We will use a common example value for kappa to demonstrate the calculation.
    kappa = 10.0
    mu = 1.0
    M = kappa

    print(f"Calculating optimal step sizes for M = {M} and mu = {mu}")
    print("=" * 50)

    # The optimal step sizes gamma_1 and gamma_2 are given by the formulas:
    # gamma_{1,2} = (4*(M+mu) +- 2*sqrt(2)*(M-mu)) / (M^2 + 6*M*mu + mu^2)

    # First, let's calculate the common denominator.
    # Denominator = M^2 + 6*M*mu + mu^2
    denominator_val = M**2 + 6 * M * mu + mu**2
    print(f"The denominator (M^2 + 6*M*mu + mu^2) is: {M**2} + {6*M*mu} + {mu**2} = {denominator_val}")

    # Next, we calculate the two parts of the numerator.
    # The first term is 4*(M+mu)
    num_term1 = 4 * (M + mu)
    print(f"The first term of the numerator (4*(M+mu)) is: 4 * ({M} + {mu}) = {num_term1}")

    # The second term is 2*sqrt(2)*(M-mu)
    num_term2 = 2 * math.sqrt(2) * (M - mu)
    print(f"The second term of the numerator (2*sqrt(2)*(M-mu)) is: 2*sqrt(2) * ({M} - {mu}) = {num_term2}")
    
    # Now, we combine them to get the numerators for gamma_1 and gamma_2.
    # Usually, gamma_1 is chosen to be the smaller step size.
    numerator_1 = num_term1 - num_term2
    numerator_2 = num_term1 + num_term2
    print(f"\nThe numerator for gamma_1 is: {num_term1} - {num_term2} = {numerator_1}")
    print(f"The numerator for gamma_2 is: {num_term1} + {num_term2} = {numerator_2}")
    print("=" * 50)

    # Finally, we calculate the values of gamma_1 and gamma_2.
    gamma_1 = numerator_1 / denominator_val
    gamma_2 = numerator_2 / denominator_val

    print(f"The optimal step size gamma_1 = {numerator_1} / {denominator_val} = {gamma_1}")
    print(f"The optimal step size gamma_2 = {numerator_2} / {denominator_val} = {gamma_2}")
    print(f"\nThe best choice for the pair (gamma_1, gamma_2) is ({gamma_1}, {gamma_2})")

if __name__ == '__main__':
    solve_optimal_steps()