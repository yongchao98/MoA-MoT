def solve_particle_problem():
    """
    Calculates contact numbers for particles in a jammed system
    and identifies the name of the unstable particles.
    """
    # The dimension of the system is given
    d = 3
    print(f"The system is in d = {d} dimensions.")
    print("-" * 40)

    # For most particles, they are jammed and have z=2d contacts
    z_jammed_factor = 2
    z_jammed = z_jammed_factor * d
    print("For a typical jammed particle, the average number of contacts is z = 2 * d.")
    print(f"Calculation: z = {z_jammed_factor} * {d} = {z_jammed}")
    print("-" * 40)

    # For the special, unstable particles, they have d+1 contacts
    z_unstable_addend = 1
    z_unstable = d + z_unstable_addend
    print("For the unstable particles, the typical number of contacts is d + 1.")
    print(f"Calculation: contacts = {d} + {z_unstable_addend} = {z_unstable}")
    print("-" * 40)

    # The name of these particles is "rattlers"
    particle_name = "rattlers"
    print(f"The name of the particles with {z_unstable} contacts that are caged but not strictly jammed is: {particle_name}.")

if __name__ == "__main__":
    solve_particle_problem()

<<<rattlers>>>