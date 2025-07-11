import sys

def solve_particle_name():
    """
    Calculates the number of contacts for jammed particles and rattlers
    in a 3D hard-sphere system and prints the name of the unstable particles.
    """
    # The dimension of the system
    d = 3

    # For most particles in a jammed state, the number of contacts is z = 2d
    z_jammed = 2 * d

    # For the special, unstable particles, the number of contacts is d + 1
    contacts_special = d + 1

    # The name for particles that are not fully jammed and have d+1 contacts
    particle_name = "rattlers"

    print(f"In a hard-sphere system in d = {d} dimensions:")
    print(f"Most jammed particles are part of the rigid structure and have z = 2 * {d} = {z_jammed} contacts.")
    print(f"The particles that are not strictly jammed and have only d + 1 = {d} + 1 = {contacts_special} contacts are called '{particle_name}'.")

solve_particle_name()