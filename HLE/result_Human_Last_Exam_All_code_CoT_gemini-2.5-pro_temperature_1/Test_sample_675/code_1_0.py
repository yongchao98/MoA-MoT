import sys

def solve_particle_name():
    """
    Calculates contact numbers for particles in a jammed system and identifies the unstable ones.
    """
    # Dimension of the system as per the problem description
    d = 3

    # The number of contacts for a particle to be mechanically stable (isostatic)
    # in a system of frictionless hard spheres is z = 2d.
    z_isostatic = 2 * d

    # The number of contacts for the unstable particles described.
    z_unstable = d + 1

    print(f"In a {d}-dimensional system of hard spheres at jamming:")
    print(f"Most particles are 'jammed' and part of the rigid backbone.")
    print(f"They satisfy the isostatic condition, having on average z = 2 * {d} = {z_isostatic} contacts.")
    print("\nHowever, a small fraction of particles are not rigidly held in place.")
    print(f"These particles have fewer contacts, typically z = {d} + 1 = {z_unstable}.")
    print("This under-coordination means they are not part of the rigid backbone and can move or 'rattle' in the cages formed by their neighbors.")
    print("\nThe specific name for these particles is 'Rattlers'.")

if __name__ == "__main__":
    solve_particle_name()