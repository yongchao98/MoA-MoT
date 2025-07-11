import numpy as np

def calculate_critical_correlation():
    """
    This function calculates the critical amount of correlation (covariance)
    required to balance potentiation and depression for selectivity to emerge.

    The plan is as follows:
    1. Define the average rate of activation, mu, as given in the problem.
    2. Apply the derived formula for the critical covariance C_sv.
    3. Print the equation, the values used, and the final result.
    """

    # The average rate of activation 'mu' is a parameter of the system.
    # As it's not specified numerically, we will use a representative value.
    # The user can change this value to explore different scenarios.
    mu = 0.1

    # Based on a stability analysis of the weight dynamics, the critical
    # condition for the emergence of selective receptive fields is determined.
    # This critical point balances the synaptic forces of potentiation and depression
    # for the selective mode of weight changes.
    # The derived formula for this critical covariance (C) is:
    # C = mu - mu^2
    critical_covariance = mu - mu**2

    # Print the final result, showing each number in the final equation.
    print("The critical amount of correlation, expressed as covariance (C), required to balance potentiation and depression for the emergence of selectivity is given by the formula:")
    print("C = mu - mu^2")
    print("\nFor the given value of mu:")
    print(f"mu = {mu}")
    print("\nThe calculation is:")
    # Using np.round to avoid potential floating point representation issues in the printout
    print(f"C = {mu} - {mu}**2")
    print(f"C = {mu} - {np.round(mu**2, 4)}")
    print(f"C = {np.round(critical_covariance, 4)}")


if __name__ == "__main__":
    calculate_critical_correlation()