def solve_particle_name():
    """
    Calculates contact numbers and identifies the particle type
    in a jammed hard-sphere system.
    """
    # Dimension of the system as given in the problem
    d = 3

    # The problem states that jammed particles have, on average, z = 2d contacts.
    # Let's calculate this value.
    z_jammed = 2 * d

    # The problem also describes special particles with z = d+1 contacts.
    # Let's calculate this value.
    z_special = d + 1

    print("In a jammed system of hard spheres in d=3 dimensions:")
    print(f"A typical jammed particle is isostatic and has an average of z = 2 * {d} = {z_jammed} contacts.")
    print(f"The special particles described are under-constrained and have z = {d} + 1 = {z_special} contacts.")
    print("This is the minimum number of contacts required to trap, or 'cage', a particle.")
    print("\nParticles that are caged but not part of the rigid, force-bearing structure are called:")
    print("Rattlers")

solve_particle_name()