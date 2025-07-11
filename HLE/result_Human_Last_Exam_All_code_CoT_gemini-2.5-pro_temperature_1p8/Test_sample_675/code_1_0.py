def solve_jamming_puzzle():
    """
    Calculates contact numbers in a 3D jammed system and identifies
    the special particles described.
    """
    # Dimension of the system
    d = 3

    # Calculate contacts for isostatic (jammed) particles
    jammed_contacts = 2 * d
    
    # Calculate contacts for the unstable particles
    unstable_contacts = d + 1
    
    # The name of these particles
    particle_name = "rattlers"

    print(f"In a d={d} dimensional hard-sphere system at jamming:")
    print(f"Most stable particles have z = 2 * {d} = {jammed_contacts} contacts.")
    print(f"The unstable particles that are not strictly jammed have z = {d} + 1 = {unstable_contacts} contacts.")
    print("\nThe name of these particles is:")
    print(particle_name)

# Execute the function
solve_jamming_puzzle()