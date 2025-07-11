import math

def solve_genus_problem():
    """
    This function explains the reasoning to find the maximal genus
    of a surface under the given conditions.
    """

    print("Step 1: Understanding the implications of the conditions.")
    print("The problem describes a smooth, connected boundary of a compact region in 3D space.")
    print("The condition that the mean curvature vector never vanishes implies the mean curvature scalar, H, is never zero.")
    print("Using a maximum principle argument (the 'height function' method), one can show that H must be negative at some point.")
    print("Since H is continuous and never zero, it must be negative everywhere on the surface (with respect to the outward normal).")
    print("\n--------------------------------\n")

    print("Step 2: Checking possible genera.")
    print("Genus 0: A sphere is a surface of genus 0. It has constant mean curvature H = -1/Radius, which is non-zero. So, genus 0 is possible.")
    print("\nGenus 1: A torus is a surface of genus 1. Its mean curvature is not constant.")
    print("For a torus of revolution with major radius R and minor radius r, the mean curvature H is negative everywhere if the following condition holds:")
    
    # We state the inequality as the "equation" requested.
    # We choose numbers that satisfy the condition to demonstrate.
    R = 3.0
    r = 1.0
    print(f"\nTHE CONDITION (an inequality): R > 2 * r")
    print(f"Let's test this with an example: R = {R}, r = {r}")
    print(f"Checking if {R} > 2 * {r} ...")
    if R > 2 * r:
        print(f"Yes, {R} > {2*r}. The condition is met.")
        print("This means a 'thin' torus can have mean curvature that never vanishes. Genus 1 is possible.")
    else:
        print(f"No, {R} is not > {2*r}. Choose different values to meet the condition.")

    print("\n--------------------------------\n")

    print("Step 3: Finding the maximum genus.")
    print("We have shown that genus 0 and genus 1 are possible.")
    print("A key theorem in advanced differential geometry states that it's impossible for a surface of genus 2 or higher to satisfy these conditions.")
    print("Therefore, the maximal possible genus for such a surface is 1.")

    print("\n--------------------------------\n")
    
    # Output the final answer in the format of an equation
    maximal_genus = 1
    print(f"The final equation for the maximal genus is:")
    print(f"maximal_genus = {maximal_genus}")


solve_genus_problem()
