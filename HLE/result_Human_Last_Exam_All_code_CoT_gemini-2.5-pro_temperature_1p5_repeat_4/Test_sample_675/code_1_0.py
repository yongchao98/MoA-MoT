import sys

def identify_particle_type(dimension):
    """
    Calculates contact numbers for particles in a jammed hard-sphere system
    and identifies the name for the unstable particles.
    """
    # Dimension of the system
    d = dimension

    # For most particles in a jammed state, the number of contacts (z)
    # follows the isostatic condition, z = 2d.
    z_jammed = 2 * d

    # For the unstable particles that are not strictly jammed,
    # the number of contacts is typically d + 1.
    z_unstable = d + 1
    
    # Print the description and calculations
    print(f"In a d={d} dimensional system:")
    print(f"Most jammed particles have z = 2 * {d} = {z_jammed} contacts.")
    print(f"The unstable particles with back-and-forth movement have z = {d} + 1 = {z_unstable} contacts.")
    print("\nThe name for these particles with an unstable number of contacts is 'Rattlers'.")

# Set the dimension as specified in the problem
d_space = 3
identify_particle_type(d_space)