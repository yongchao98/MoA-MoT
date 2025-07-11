import sys

def identify_particles():
    """
    This script calculates the number of contacts for jammed and unstable
    particles in a 3D hard-sphere system and identifies the unstable ones.
    """
    # Dimension of the system as provided in the problem
    d = 3

    # For most particles, the system is isostatic (rigidly jammed).
    # The number of contacts, z, is 2*d.
    isostatic_contacts = 2 * d

    # For the unstable particles, the number of contacts is d+1.
    unstable_contacts = d + 1
    
    # The name of these particles that are not part of the rigid backbone
    # and can move within their cages.
    particle_name = "Rattlers"

    print("In a hard-sphere system at jamming in dimension d=3:")
    print(f"Most particles are isostatic and have z = 2 * {d} = {isostatic_contacts} contacts.")
    print(f"The particles with an unstable contact have z = {d} + 1 = {unstable_contacts} contacts.")
    print("\nThe name for these unstable particles is:")
    print(particle_name)

if __name__ == "__main__":
    identify_particles()