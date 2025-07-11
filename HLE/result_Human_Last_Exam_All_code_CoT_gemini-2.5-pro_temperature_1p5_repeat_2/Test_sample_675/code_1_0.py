import sys

def find_particle_name():
    """
    Calculates particle contacts in a jammed hard-sphere system
    and identifies the name of the unstable particles.
    """
    # Dimension of the system
    d = 3

    # On average, each jammed particle has z = 2d contacts
    jammed_contacts = 2 * d

    # The special, unstable particles have d+1 contacts
    unstable_contacts = d + 1

    particle_name = "rattlers"

    # Print the final explanation including the calculated numbers and the name
    print(f"In a {d}-dimensional system, a jammed particle has {jammed_contacts} contacts on average ({d}*2).")
    print(f"The special particles, which are not strictly jammed and can move back-and-forth, have {d} + 1 = {unstable_contacts} contacts.")
    print(f"The name for these particles is: {particle_name}")

if __name__ == "__main__":
    find_particle_name()
