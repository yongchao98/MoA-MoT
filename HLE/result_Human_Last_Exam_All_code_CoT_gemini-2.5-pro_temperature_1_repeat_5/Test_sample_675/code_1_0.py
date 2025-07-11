import sys

def solve():
    """
    This script calculates the number of contacts for jammed and special particles
    in a 3D hard-sphere system and identifies the name of the special particles.
    """
    # The dimension of the system
    d = 3

    # Calculate the average number of contacts for a jammed particle
    # The equation is z = 2d
    z_jammed = 2 * d

    # Calculate the number of contacts for the special, unstable particle
    # The equation is z = d + 1
    z_special = d + 1

    print(f"In a system of dimension d = {d}:")
    print(f"The average number of contacts for a jammed particle is 2 * {d} = {z_jammed}.")
    print(f"The special, unstable particles have {d} + 1 = {z_special} contacts.")

    # Identify the name based on the description
    particle_name = "Rattlers"
    print(f"\nParticles that are caged but not strictly jammed, and exhibit a back-and-forth movement, are called: {particle_name}")

    # Output the final answer in the specified format
    # The double braces {{}} are used to escape the f-string formatting for the outer braces
    # to literally print <<<Rattlers>>>
    sys.stdout.write(f"\n<<<{particle_name}>>>\n")

solve()