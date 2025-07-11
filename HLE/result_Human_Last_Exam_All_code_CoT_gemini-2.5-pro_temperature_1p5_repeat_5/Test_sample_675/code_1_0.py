import sys

def find_particle_name():
    """
    Calculates contact numbers for particles in a jammed system and identifies the name of the unstable particles.
    """
    # The dimension of the system as given in the problem
    d = 3

    # Calculate the average number of contacts for a stable, jammed particle (isostatic condition)
    # The equation is z = 2 * d
    z_stable = 2 * d

    # Calculate the number of contacts for the unstable particles
    # The equation is contacts = d + 1
    contacts_unstable = d + 1

    particle_name = "rattlers"

    print("Analyzing the properties of particles in a hard sphere system:")
    print(f"Given dimension d = {d}")
    
    print("\nFor a stable, jammed particle, the number of contacts 'z' is given by the equation: z = 2 * d")
    print(f"Calculation: z = 2 * {d} = {z_stable}")
    
    print("\nFor the unstable particle described, the number of contacts is given by the equation: contacts = d + 1")
    print(f"Calculation: contacts = {d} + 1 = {contacts_unstable}")
    
    print(f"\nParticles that are caged by their neighbors but are not strictly jammed and have {contacts_unstable} contacts (fewer than the stable {z_stable} contacts) are called '{particle_name}'.")
    print("Their unstable nature allows them to 'rattle' within their cages.")


if __name__ == "__main__":
    find_particle_name()
