def solve_particle_name():
    """
    Calculates contact numbers for particles in a 3D hard sphere system
    and identifies the name of the unstable particles.
    """
    # Dimension of the system
    d = 3
    print(f"In a system of dimension d = {d}:")

    # Calculate average contacts for a jammed particle (isostatic condition)
    jammed_contacts_factor = 2
    jammed_contacts = jammed_contacts_factor * d
    print(f"The average number of contacts for a stable, jammed particle is z = {jammed_contacts_factor} * {d} = {jammed_contacts}.")

    # Calculate contacts for the unstable particle
    unstable_contacts_addend = 1
    unstable_contacts = d + unstable_contacts_addend
    print(f"The number of contacts for the unstable particles is d + 1 = {d} + {unstable_contacts_addend} = {unstable_contacts}.")

    # Define the name of these particles
    particle_name = "rattlers"
    print(f"\nParticles with {unstable_contacts} contacts that are not strictly jammed and can move are called: {particle_name}.")

solve_particle_name()