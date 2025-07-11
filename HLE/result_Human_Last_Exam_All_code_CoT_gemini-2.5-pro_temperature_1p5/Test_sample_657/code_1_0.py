def print_stabilizing_controller():
    """
    This function prints the components of the transfer function for the set of
    all proper stabilizing controllers H_2(s) for the given system.

    The controller H_2(s) is parametrized by a function K(s),
    which must be stable and proper.
    """

    # The numerator of H_2(s) as a function of s and K(s)
    # Numerator = (4*s**2 + 8*s + 4) + (s**2 - 1)*K(s)
    numerator_str = "(4*s**2 + 8*s + 4) + (s**2 - 1)*K(s)"

    # The denominator of H_2(s) as a function of s and K(s)
    # Denominator = (s**2 - 1) - s*K(s)
    denominator_str = "(s**2 - 1) - s*K(s)"

    print("The set of all proper stabilizing controllers H_2(s) is given by:")
    print("H_2(s) = Numerator / Denominator\n")
    print(f"Numerator: {numerator_str}")
    print(f"Denominator: {denominator_str}\n")
    print("Where K(s) is any arbitrary stable and proper transfer function.")

# Execute the function to print the result
print_stabilizing_controller()
