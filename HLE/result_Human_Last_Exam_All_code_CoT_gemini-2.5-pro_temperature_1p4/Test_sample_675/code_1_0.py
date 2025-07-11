def find_particle_name():
    """
    Calculates the number of contacts for jammed and unstable particles
    in a 3D hard-sphere system and identifies the unstable particles.
    """
    # The dimension of the system
    d = 3

    # For most jammed particles, the average number of contacts is z = 2d
    jammed_contacts = 2 * d

    # For the unstable particles, the number of contacts is d + 1
    unstable_contacts = d + 1

    print("In a hard-sphere system at jamming:")
    print(f"The dimension (d) is: {d}")
    print(f"Most jammed particles have z = 2 * d = 2 * {d} = {jammed_contacts} contacts.")
    print(f"The special, unstable particles have d + 1 = {d} + 1 = {unstable_contacts} contacts.")
    print("\nThese under-coordinated particles, which are not strictly jammed and can move back-and-forth, are called:")
    print("rattlers")

# Execute the function to find the answer
find_particle_name()
<<<rattlers>>>