def solve_particle_name():
    """
    Calculates the number of contacts for jammed particles and rattlers
    in a 3D hard-sphere system and identifies the name of the latter.
    """
    # The dimension of the system
    d = 3

    # For most jammed particles, the number of contacts is z = 2d
    z_jammed = 2 * d

    # For the specific particles in question (rattlers), the number of contacts is z = d + 1
    z_rattler = d + 1

    # The name of these particles
    particle_name = "rattlers"

    print(f"In a system with dimension d = {d}:")
    print(f"The average number of contacts for a jammed particle is 2 * {d} = {z_jammed}.")
    print(f"The number of contacts for the unstable particles is {d} + 1 = {z_rattler}.")
    print(f"These particles are called: {particle_name}")

solve_particle_name()

<<<rattlers>>>