import numpy as np

def calculate_critical_correlation(N_v, N_s, mu, theta):
    """
    Calculates the critical amount of correlation (C_vs) required to balance
    potentiation and depression in the network.

    Args:
        N_v (int): Number of neurons in the input layer v.
        N_s (int): Number of neurons in the input layer s.
        mu (float): The average rate of activation for input neurons.
        theta (float): The heterosynaptic offset constant.

    Returns:
        float: The critical correlation C_vs = <s_k * v_k>.
    """
    
    # The derived formula for the critical correlation
    C_vs = (N_v + N_s) * mu * (theta - mu) - mu * (1 - 2 * mu)
    
    return C_vs

def main():
    """
    Main function to demonstrate the calculation with example values.
    """
    # Example parameters
    # Number of neurons in each input layer
    N_v = 100
    N_s = 100
    # Average activation rate (probability of firing in a small time bin)
    mu = 0.01
    # Heterosynaptic depression threshold
    theta = 0.015

    print("This script calculates the critical correlation 'C_vs' required to balance synaptic potentiation and depression.")
    print("The formula is derived from the stability condition <r_i * v_k> = theta * <r_i>.\n")

    print("Formula:")
    print("C_vs = (N_v + N_s) * mu * (theta - mu) - mu * (1 - 2*mu)\n")

    print("Using the following example values:")
    print(f"N_v = {N_v}")
    print(f"N_s = {N_s}")
    print(f"mu = {mu}")
    print(f"theta = {theta}\n")

    # Calculate the critical correlation
    critical_correlation = calculate_critical_correlation(N_v, N_s, mu, theta)

    print("Plugging the numbers into the equation:")
    # Using f-string to display the equation with the numbers
    equation_str = f"C_vs = ({N_v} + {N_s}) * {mu} * ({theta} - {mu}) - {mu} * (1 - 2 * {mu})"
    print(equation_str)
    
    # Calculate intermediate steps for clarity
    term1 = (N_v + N_s) * mu * (theta - mu)
    term2 = mu * (1 - 2 * mu)
    print(f"C_vs = {term1:.6f} - {term2:.6f}")
    
    print(f"\nThe critical amount of correlation, C_vs, is: {critical_correlation:.6f}")

    # For context, compare to the correlation of independent neurons
    independent_corr = mu**2
    print(f"For comparison, the correlation for independent neurons (mu^2) would be: {independent_corr:.6f}")


if __name__ == "__main__":
    main()
    print("\nFinal Answer Formula:")
    print("C_vs = (N_v + N_s)*μ*(θ - μ) - μ*(1 - 2*μ)")
    # The user is asked for 'the' critical amount, so we return the symbolic formula
    # rather than a specific number, as the numbers are just examples.
    # However, to be helpful to the user and in case they wanted a numerical answer from the example,
    # the code above calculates it.
    # Following the thought process, the critical correlation is the condition for balance.
    # We output the formula as the most general answer.
    # Wrapping final symbolic answer in <<<>>>
    final_answer_formula = "(N_v + N_s)*mu*(theta - mu) - mu*(1 - 2*mu)"


<<< (N_v + N_s)*μ*(θ - μ) - μ*(1 - 2*μ) >>>