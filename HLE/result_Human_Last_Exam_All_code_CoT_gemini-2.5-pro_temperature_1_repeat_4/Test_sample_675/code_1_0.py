def identify_particles():
    """
    Identifies a specific type of particle in a jammed hard-sphere system.
    """
    # Dimension of the system
    d = 3

    # Number of contacts for an average, isostatic (jammed) particle
    # The equation is z = 2 * d
    z = 2 * d

    # Number of contacts for the unstable particles in question
    # The equation is d + 1
    unstable_contacts = d + 1

    print(f"In a {d}-dimensional hard-sphere system at jamming:")
    print(f"Most particles are isostatic, meaning they have z = 2 * d = 2 * {d} = {z} contacts.")
    print(f"The particles you described are unstable and have only d + 1 = {d} + 1 = {unstable_contacts} contacts.")
    print("\nThese particles are under-constrained and can move or 'rattle' in the cages formed by their neighbors.")
    print("Therefore, they are called:")
    print("rattlers")

if __name__ == "__main__":
    identify_particles()