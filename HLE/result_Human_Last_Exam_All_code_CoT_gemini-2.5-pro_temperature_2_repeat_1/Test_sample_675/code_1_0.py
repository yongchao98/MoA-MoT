def solve_particle_name():
    """
    Calculates the number of contacts for isostatic and rattler particles
    in a 3D hard-sphere system and identifies the name of the rattlers.
    """
    # The dimension of the system
    d = 3

    # For most particles at jamming, the number of contacts is z = 2d
    # This is the isostatic condition for mechanical stability.
    isostatic_contacts = 2 * d

    # For the special, unstable particles described, the number of contacts is d+1.
    # These particles are under-constrained.
    unstable_contacts = d + 1

    print(f"In a d={d} system, a typically jammed particle has z = 2 * {d} = {isostatic_contacts} contacts.")
    print(f"The particles in question, having unstable contacts and moving back-and-forth, have z = {d} + 1 = {unstable_contacts} contacts.")
    print("These particles are not part of the main rigid structure and are called:")
    print("Rattlers")

solve_particle_name()