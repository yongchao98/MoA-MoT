def identify_particle_type():
    """
    Identifies the name of special particles in a jammed hard-sphere system
    based on their number of contacts.
    """
    # Dimension of the system
    d = 3

    # The formula for the average number of contacts for a typical, isostatically jammed particle
    # z_jammed = 2 * d
    num_contacts_jammed = 2 * d

    # The formula for the number of contacts for the unstable particles
    # z_special = d + 1
    num_contacts_special = d + 1

    print(f"In a {d}-dimensional hard-sphere system at jamming:")
    print(f"A typical, fully jammed particle has an average of 2*d = 2*{d} = {num_contacts_jammed} contacts.")
    print(f"The special, unstable particles described have d+1 = {d}+1 = {num_contacts_special} contacts.")
    print("\nThese particles are not part of the rigid backbone of the material and can move locally.")
    print("The name of these particles is: Rattlers")

# Run the function to display the answer
identify_particle_type()