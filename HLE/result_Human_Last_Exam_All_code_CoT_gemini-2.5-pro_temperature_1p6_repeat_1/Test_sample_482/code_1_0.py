import math

def calculate_critical_correlation(mu, theta):
    """
    Calculates the critical amount of correlation required to balance potentiation and depression.

    Args:
        mu (float): The average rate of activation for input neurons.
        theta (float): The heterosynaptic offset constant.

    Returns:
        float: The critical correlation C = E[v*s], or None if parameters are invalid.
    """
    # The derivation C = 2*mu*theta - mu is physically meaningful only if C >= 0,
    # which implies theta >= 0.5. We'll also check if mu and theta are in a valid range.
    if not (0 < mu < 1):
        print("Error: The average rate 'mu' must be between 0 and 1.")
        return None
    if not (0 < theta < 1):
        print("Error: The threshold 'theta' must be between 0 and 1.")
        return None
    if theta < 0.5:
        print(f"Warning: With theta={theta} < 0.5, a cooperative stable state is not possible.")
        print("The calculated correlation would be negative, which is not physically meaningful.")
        # We calculate it anyway as per the formula.
    
    # The critical amount of correlation C = E[v*s]
    critical_correlation = 2 * mu * theta - mu
    
    print("The critical amount of correlation 'C' is determined by the equation: C = 2 * mu * theta - mu")
    print("This condition allows for a stable, non-trivial wiring structure to form in the network.")
    print("\nUsing the provided values:")
    print(f"mu (average activation rate) = {mu}")
    print(f"theta (heterosynaptic threshold) = {theta}")
    print("\nThe final equation with numbers is:")
    print(f"{critical_correlation:.4f} = 2 * {mu} * {theta} - {mu}")
    
    return critical_correlation

# --- Parameters ---
# mu: Average rate of activation for inputs v and s. Must be in (0, 1).
# We choose a value representative of sparse activity.
mu = 0.1

# theta: Heterosynaptic offset constant. Must be in (0, 1).
# For a cooperative solution to exist, we need theta >= 0.5.
theta = 0.6

# --- Calculation ---
C_critical = calculate_critical_correlation(mu, theta)

if C_critical is not None:
    # Final answer block as requested
    print(f"\n<<<The critical correlation is {C_critical:.4f}>>>")
