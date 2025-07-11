import numpy as np

def calculate_critical_correlation(theta, mu, V=None):
    """
    Calculates the critical correlation required to balance potentiation and depression.

    The formula is derived from the stability analysis of the weight dynamics:
    C_crit = 2 * theta * mu - V
    where:
    - C_crit: Critical correlation <s_k * v_k>
    - theta:  The heterosynaptic offset constant.
    - mu:     The average rate of activation for inputs v and s (<v_k>).
    - V:      The second moment of the input activity (<v_k^2>).

    For a Poisson process modeled as binary events (0 or 1), V is equal to mu.
    This function uses this simplification if V is not provided.
    """
    if V is None:
        # For a Poisson process modeled by a Binomial variable (event or no event),
        # where the probability of an event is mu, the second moment <v^2> is also mu.
        V = mu
        print("Parameter V (second moment of input) not provided.")
        print(f"Assuming binary input activity based on a Poisson process, so V = mu = {mu:.3f}\n")


    # Calculate the critical correlation
    c_critical = 2 * theta * mu - V

    # --- Output the results ---
    print("The formula for the critical amount of correlation (C_crit) is:")
    print("C_crit = 2 * theta * mu - V\n")

    print("Substituting the given values into the equation:")
    # Using f-string formatting to display the equation with numbers
    # Each number is explicitly shown as requested.
    print(f"C_crit = 2 * {theta:.3f} * {mu:.3f} - {V:.3f}")

    # Calculate and print the final result
    term1 = 2 * theta
    term2 = term1 * mu
    print(f"C_crit = {term1:.3f} * {mu:.3f} - {V:.3f}")
    print(f"C_crit = {term2:.3f} - {V:.3f}")
    
    # Final answer wrapped in <<<>>>
    print(f"<<<{c_critical:.3f}>>>")


if __name__ == '__main__':
    # --- Define model parameters ---
    # These are example values. You can replace them with your own.

    # heterosynaptic offset constant
    theta = 0.8

    # The problem states an inter-event interval of 150s for the Poisson process.
    # The average rate 'mu' in the learning rule depends on the time step 'dt' used
    # in the simulation (mu = dt / 150s).
    # Let's assume a simulation where mu = 0.05.
    mu = 0.05

    # We do not need to provide V as the function will set V = mu
    # for this type of input.
    calculate_critical_correlation(theta, mu)
