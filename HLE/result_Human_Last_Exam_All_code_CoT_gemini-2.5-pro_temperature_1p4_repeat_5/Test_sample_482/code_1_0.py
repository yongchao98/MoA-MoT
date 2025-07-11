import math

def calculate_critical_correlation(N, mu, theta):
    """
    Calculates the critical amount of correlation between input populations v and s.

    This value of correlation ensures that the potentiation and depression effects
    in the synaptic weights are balanced, allowing for a stable and symmetric
    equilibrium where the total synaptic strengths from both populations are equal.

    The formula is derived from the stability analysis of the weight dynamics.
    C_sv = 2*theta*mu - mu/N - (1 - 1/N)*mu^2

    Args:
        N (int): The number of neurons in each input population (Nv = Ns = N).
        mu (float): The average rate of activation for any input neuron.
        theta (float): The heterosynaptic offset constant from the learning rule.

    Returns:
        float: The critical correlation C_sv.
    """
    # Check if parameters lead to a positive correlation, which is required
    # for the correlation to be physically meaningful (as a joint probability).
    # This requires C_sv > 0, which approximately means 2*theta > mu for large N.
    if mu * (2 * theta - mu) < 0:
        print("Warning: Parameters may result in a negative correlation, "
              "which is not physically meaningful. Check if 2*theta > mu.")

    # Calculate the critical correlation
    c_sv = 2 * theta * mu - mu / N - (1 - 1/N) * (mu**2)

    # Print the equation with the given numbers
    print("Calculating the critical correlation C_sv:")
    print(f"C_sv = 2 * {theta} * {mu} - {mu} / {N} - (1 - 1/{N}) * {mu}^2")
    print("-" * 20)
    print(f"Final calculated value:")
    print(f"C_sv = {c_sv}")
    
    return c_sv

if __name__ == '__main__':
    # --- Parameters ---
    # Since specific values were not provided in the problem description,
    # we use the following example values.
    # N: Number of neurons in each input layer.
    N = 100
    # mu: Average activation rate of an input neuron (e.g., probability of firing in a time bin).
    mu = 0.01
    # theta: Heterosynaptic offset constant. Should be of a similar order of magnitude as mu.
    theta = 0.1

    print(f"Given parameters:\n  N = {N}\n  mu = {mu}\n  theta = {theta}\n")

    # Calculate and print the result
    critical_correlation = calculate_critical_correlation(N, mu, theta)
    # The final answer format as requested.
    print(f"\n<<< {critical_correlation} >>>")