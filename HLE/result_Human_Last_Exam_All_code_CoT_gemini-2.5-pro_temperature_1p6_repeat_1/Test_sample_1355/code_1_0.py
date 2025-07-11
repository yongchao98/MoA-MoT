def calculate_conductance_moment_ratio():
    """
    Calculates and explains the ratio of the fourth to the first statistical
    moment of conductance for a disordered Majorana wire at the critical point.
    """
    
    # In the theory of disordered topological systems, an ensemble of Majorana wires
    # at the critical point between the trivial and topological phases belongs to the
    # BDI symmetry class. In the limit of long wire length, the statistical
    # properties of the dimensionless conductance, g, become universal.
    
    # A key result for this universality class is that the n-th statistical moment
    # of the conductance, <g^n>, is given by a remarkably simple formula:
    # <g^n> = 1 / (n + 1).
    # This implies that the conductance effectively behaves as a random variable
    # uniformly distributed on the interval [0, 1].

    # We need the fourth statistical moment (n=4):
    n_fourth = 4
    g_fourth_moment = 1 / (n_fourth + 1)

    # We also need the average value, which is the first statistical moment (n=1):
    n_first = 1
    g_average = 1 / (n_first + 1)

    # The problem asks for the ratio between these two quantities.
    ratio = g_fourth_moment / g_average

    # Print the explanation and the final equation with numerical values as requested.
    print("For a disordered Majorana wire at the critical point, the statistical moments of the dimensionless conductance 'g' are given by <g^n> = 1/(n+1).")
    
    print("\nThe average value (first moment, n=1) of the conductance is:")
    print(f"<g> = 1 / ({n_first} + 1) = {g_average}")
    
    print("\nThe fourth statistical moment (n=4) of the conductance is:")
    print(f"<g^4> = 1 / ({n_fourth} + 1) = {g_fourth_moment}")
    
    print("\nThe ratio between the fourth moment and the average value is therefore:")
    print(f"<g^4> / <g> = {g_fourth_moment} / {g_average} = {ratio}")

if __name__ == "__main__":
    calculate_conductance_moment_ratio()