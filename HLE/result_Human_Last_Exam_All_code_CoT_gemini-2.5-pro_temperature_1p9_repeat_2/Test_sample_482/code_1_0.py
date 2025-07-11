import argparse

def solve_critical_correlation(mu):
    """
    Calculates the 'critical amount of correlation' C required to balance
    potentiation and depression in the described neural network.

    Args:
        mu (float): The average rate of activation of the input neurons (a value
                    between 0 and 1).
    """
    # The 'critical amount of correlation' is the covariance C between corresponding
    # input neurons (e.g., v_k and s_k) that allows the network to develop
    # stable, selective weights. For this to happen, C must balance the synaptic
    # updates. The analytical derivation shows that this occurs when the covariance C
    # is equal to the variance of the presynaptic input activity.

    # Check if the input mu is a valid probability
    if not (0 <= mu <= 1):
        print("Error: The average rate 'mu' must be a probability, i.e., between 0 and 1.")
        return

    # For an input modeled as a Bernoulli random variable with probability mu,
    # the variance is given by the formula: mu * (1 - mu).
    critical_correlation = mu * (1 - mu)

    # Output the explanation and the final equation with the numbers plugged in.
    print("The analysis reveals that the critical amount of correlation, C, is equal to")
    print("the variance of the presynaptic neuron's activity.")
    print("\nThe governing equation for the critical correlation (C) is:")
    print("C = mu * (1 - mu)")
    print("-" * 50)
    print(f"Given an average input rate mu = {mu}:")
    
    # Final equation with the provided numbers, as requested.
    print("\nThe final calculation is:")
    print(f"C = {mu} * (1 - {mu}) = {critical_correlation:.4f}")
    print("-" * 50)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Determine the critical amount of correlation C in a feedforward neural network.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "--mu",
        type=float,
        default=0.1,
        help="The average rate of activation (mu) for the input neurons.\n"
             "It represents the probability of firing in a time step.\n"
             "Default value is 0.1."
    )
    args = parser.parse_args()
    solve_critical_correlation(args.mu)