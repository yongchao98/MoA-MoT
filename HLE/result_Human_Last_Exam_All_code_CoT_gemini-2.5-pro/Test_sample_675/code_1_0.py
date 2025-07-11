def calculate_contacts():
    """
    Calculates and prints the number of contacts for jammed and unstable
    particles in a 3D hard-sphere system.
    """
    # Dimension of the system
    d = 3

    # Calculate the average number of contacts for a typically jammed particle
    # The equation is z = 2d
    jammed_contacts = 2 * d

    # Calculate the number of contacts for the unstable particles
    # The equation is z = d + 1
    unstable_contacts = d + 1

    # Print the results, showing the numbers used in the equations
    print(f"In a d={d} system, most jammed particles have z = 2 * {d} = {jammed_contacts} contacts.")
    print(f"The special, unstable particles have z = {d} + 1 = {unstable_contacts} contacts.")

calculate_contacts()