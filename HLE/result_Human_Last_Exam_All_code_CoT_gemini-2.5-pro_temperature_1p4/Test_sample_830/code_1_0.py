import math

def calculate_limit_value(q, g, h):
    """
    Calculates the value of the limit based on the derived formula.

    The problem is to find the value of:
    liminf_{N(I) -> +inf} (phi(I) * log(log(N(I)))) / N(I)

    The derived value is (e^-gamma * q^g * log(q)) / h.

    Args:
        q (int): The characteristic of the local field, must be > 0.
        g (int): The genus of the function field.
        h (int): The ideal class number of the ring R.
    """
    # The Euler-Mascheroni constant
    gamma = 0.57721566490153286060651209008240243104215933593992

    if q <= 1:
        print("Error: q must be a prime power greater than 1.")
        return
    if h <= 0:
        print("Error: h (class number) must be positive.")
        return

    # Calculate the components of the formula
    e_to_neg_gamma = math.exp(-gamma)
    q_to_g = q**g
    log_q = math.log(q)

    # The final value of the limit
    result = (e_to_neg_gamma * q_to_g * log_q) / h
    
    # Print the equation with the values used
    print(f"For q = {q}, g = {g}, h = {h}:")
    print("The formula for the limit is: (e^-gamma * q^g * log(q)) / h")
    print(f"e^-gamma = {e_to_neg_gamma}")
    print(f"q^g = {q_to_g}")
    print(f"log(q) = {log_q}")
    print(f"h = {h}")
    print("\nFinal Equation:")
    print(f"({e_to_neg_gamma} * {q_to_g} * {log_q}) / {h} = {result}")

# Example usage with some placeholder values
# You can change these to any valid values for q, g, and h.
q_val = 5
g_val = 2
h_val = 10

calculate_limit_value(q_val, g_val, h_val)