import sympy as sp

def print_critical_correlation_equation():
    """
    This function defines the terms for input statistics symbolically
    and prints the derived equation for the critical correlation.
    """

    # Define the symbols for the statistical properties of the neural activities.
    # C_vs represents the covariance between corresponding neurons in input populations v and s.
    Cov_vs = sp.Symbol("Cov(v,s)")

    # Var_v represents the variance of a neuron in input population v.
    # We assume v and s have identical statistics, so Var(v) = Var(s).
    Var_v = sp.Symbol("Var(v)")

    # The theoretical analysis shows that the balance between potentiation and depression
    # in the network leads to a critical point. This critical point is defined by a specific
    # relationship between the correlation of the inputs and their inherent variability.
    # The derived condition for this critical correlation is:
    critical_condition = sp.Eq(Cov_vs, Var_v)

    # Print the final equation that represents this result.
    print("The 'critical amount of correlation' is the condition required to balance potentiation and depression.")
    print("This balance is achieved when the following relationship holds:")

    # Print each part of the final equation.
    print(f"\n{critical_condition.lhs} = {critical_condition.rhs}\n")

    print("This equation states that the system reaches a critical point when the covariance between the two")
    print("input populations (v and s) becomes equal to the variance of an input population (v).")

print_critical_correlation_equation()