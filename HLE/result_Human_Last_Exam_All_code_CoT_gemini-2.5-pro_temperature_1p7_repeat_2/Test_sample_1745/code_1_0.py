def display_potential_distribution():
    """
    This function prints the final derived expression for the electrical
    double-layer potential distribution psi(y) in the microchannel.
    
    The final expression is:
    ψ(y) = z_1 * (1 + β*k) * sinh(k*(H - y)) / sinh(k*H)
    """

    print("The final expression for the Electrical double-layer potential distribution ψ(y) is:")
    
    # We print the equation component by component to fulfill the requirement
    # of outputting each number in the final equation, which is '1' in this case.
    print("ψ(y) = z_1 * (", 1, " + β*k) * sinh(k*(H - y)) / sinh(k*H)")

# Execute the function to display the result
display_potential_distribution()