def identify_particle_type():
    """
    Calculates the number of contacts for jammed and unstable particles
    in a 3D hard-sphere system and identifies the name of the unstable particles.
    """
    # Dimension of the system
    d = 3

    # For most particles, the system is isostatic (jammed).
    # The average number of contacts 'z' is 2 * d.
    z_jammed = 2 * d

    # For the unstable particles in question, the number of contacts is d + 1.
    z_unstable = d + 1

    # The name for particles that are not part of the rigid backbone of the
    # jammed solid and can move within their local cages.
    particle_name = "rattlers"

    print(f"In a {d}-dimensional system:")
    print(f"A fully jammed particle has, on average, 2 * {d} = {z_jammed} contacts.")
    print(f"The unstable particle described has {d} + 1 = {z_unstable} contacts.")
    print(f"Because these particles are under-constrained and can 'rattle' in their cages, they are called: {particle_name}.")

identify_particle_type()