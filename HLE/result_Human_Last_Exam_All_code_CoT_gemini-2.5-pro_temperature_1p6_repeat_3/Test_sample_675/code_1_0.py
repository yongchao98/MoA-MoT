def solve_particle_identity():
    """
    Calculates contact numbers and identifies the special particles in a jammed system.
    """
    # Dimension of the system
    d = 3

    # The average number of contacts for a stable, jammed particle follows the
    # isostatic condition, z = 2d.
    jammed_particle_contacts = 2 * d

    # The number of contacts for the unstable particles described in the problem.
    unstable_particle_contacts = d + 1

    # The name for these particles that are caged but not jammed.
    particle_name = "rattlers"

    print(f"In a {d}-dimensional system, a stable particle has z = 2*d = {jammed_particle_contacts} contacts on average.")
    print(f"The particles with only d+1 = {unstable_particle_contacts} contacts are not strictly jammed and are called '{particle_name}'.")

solve_particle_identity()