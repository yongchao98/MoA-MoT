import numpy as np

def calculate_critical_correlation():
    """
    Calculates the critical amount of correlation between input populations
    that balances potentiation and depression.

    The derivation starts from the weight update rule, dW/dt ∝ r(x - θ), where 'r' is the
    postsynaptic rate, 'x' is the presynaptic rate, and 'θ' is the depression threshold.
    The balance point is where the expected weight change is zero: <r(x - θ)> = 0.

    Assuming symmetric pathways for inputs v and s, the postsynaptic rate 'r' is proportional
    to the sum of the inputs: r ∝ (v + s).
    The balance condition for the v-pathway becomes <(v+s)(v-θ)> = 0.
    Expanding this gives: <v^2> - θ<v> + <sv> - θ<s> = 0.

    We model neuron activity as a binary process (0 or 1) with an average rate μ.
    This gives <v> = <s> = μ, and <v^2> = μ.
    The correlation term is C = <sv>.

    Substituting these into the equation:
    μ - θμ + C - θμ = 0
    μ + C = 2θμ
    C = 2θμ - μ
    C = μ(2θ - 1)

    This is the critical correlation C = <sv> required for balance.
    Below, we calculate this value using example parameters for μ and θ.
    """
    # --- Parameters ---
    # mu: The average firing rate of input neurons. Must be > 0.
    # Let's assume a plausible value, e.g., 0.1 (10% activity).
    mu = 0.1

    # theta: The heterosynaptic offset constant from the learning rule.
    # This value determines the balance between Hebbian potentiation and
    # heterosynaptic depression. Let's assume a value, e.g., 0.6.
    # For the correlation C to be non-negative (as <sv> >= 0), we need θ >= 0.5.
    theta = 0.6

    # --- Calculation ---
    # Using the derived formula C_crit = mu * (2 * theta - 1)
    c_crit = mu * (2 * theta - 1)

    # --- Output ---
    print("The critical amount of correlation 'C' is determined by the balance between potentiation and depression.")
    print("The derived formula for this critical correlation is: C = mu * (2 * theta - 1)")
    print("\nAssuming the following parameter values:")
    print(f"Average input rate (mu) = {mu}")
    print(f"Heterosynaptic offset (theta) = {theta}")
    print("\nWe can calculate the value of C:")
    print(f"C = {mu} * (2 * {theta} - 1)")
    print(f"C = {mu} * ({2 * theta} - 1)")
    print(f"C = {mu} * ({round(2 * theta - 1, 4)})")
    print(f"C = {round(c_crit, 4)}")
    print("\nThis value represents the required average product of activities, C = <s*v>, between the two input populations to stabilize the network weights.")

calculate_critical_correlation()
<<<C = mu * (2 * theta - 1)>>>