import math

def calculate_critical_correlation(Nv, Ns, mu, theta):
    """
    Calculates the critical amount of correlation 'C' required to balance
    potentiation and depression in the described neural network.

    Args:
        Nv (int): Number of neurons in the input layer v.
        Ns (int): Number of neurons in the input layer s.
        mu (float): The average rate of activation for v and s.
        theta (float): The heterosynaptic offset constant.
    """
    print(f"Parameters: Nv = {Nv}, Ns = {Ns}, mu = {mu}, theta = {theta}\n")

    # The derivation shows that a unique critical correlation C is found
    # when the two input populations have the same size.
    if Nv != Ns:
        print("A stable selective state is only possible under specific constraints when Nv != Ns.")
        print("This requires the correlation C to be equal to mu^2 (uncorrelated inputs),")
        print("and for theta to satisfy a separate constraint.")
        C = mu**2
        print(f"\nFor this case, the required correlation C = mu^2")
        print(f"C = {mu}**2")
        print(f"C = {C}")
        return C
    
    # When Nv = Ns = N, we can calculate the critical correlation directly.
    N = Nv
    print("Derivation for Nv = Ns = N:")
    print("C = (2 * N * theta * mu - mu - (N - 1) * mu**2) / N")
    
    # Showing the calculation step-by-step with the provided numbers
    term1 = 2 * N * theta * mu
    term2 = mu
    term3 = (N - 1) * mu**2
    numerator = term1 - term2 - term3
    C = numerator / N

    print(f"C = (2 * {N} * {theta} * {mu} - {mu} - ({N} - 1) * {mu}**2) / {N}")
    print(f"C = ({term1} - {term2} - ({N-1}) * {mu**2}) / {N}")
    print(f"C = ({term1} - {term2} - {term3}) / {N}")
    print(f"C = ({numerator}) / {N}")
    print(f"Final calculated critical correlation C = {C}")
    
    # Check if the result is physically possible (C must be <= mu)
    if C > mu:
        print("\nWarning: The calculated C is greater than mu.")
        print("This indicates that the depression term (theta) is too strong for the network to stabilize, even with maximum possible correlation.")
    elif C < 0:
        print("\nWarning: The calculated C is negative.")
        print("This indicates a requirement for anti-correlation that might be difficult to achieve biologically.")

    return C

if __name__ == '__main__':
    # --- User-configurable parameters ---
    # Number of neurons in the input layers
    N_v = 100
    N_s = 100
    # Average rate of activation (e.g., probability of firing in a time step)
    mu_val = 0.1
    # Heterosynaptic offset constant
    theta_val = 0.5
    # --- End of parameters ---

    critical_c = calculate_critical_correlation(N_v, N_s, mu_val, theta_val)
    # The final answer tag requires a single numerical value.
    # We will output the value calculated with the example parameters.
    # To use in the final answer tag, format it as a simple number.
    # For example: <<<0.0891>>>