import sys

# Define a function to solve the problem and print the steps.
def solve_surface_genus():
    """
    Calculates the smallest possible genus for a closed surface containing
    a given open surface based on topological principles.
    """
    # The genus of the initial surface, Sigma.
    g_sigma = 10

    print("Step 1: Identify the properties of the given surface Sigma.")
    print(f"The genus of Sigma is g_sigma = {g_sigma}.")
    print("Sigma has a single boundary component, which is an unknotted circle.")
    print("-" * 20)

    print("Step 2: Understand the goal.")
    print("We need to find the smallest genus 'g' for a closed surface Sigma'")
    print("such that Sigma is a subset of Sigma'.")
    print("This is equivalent to 'capping' the boundary of Sigma.")
    print("-" * 20)

    print("Step 3: Relate the genera of the surfaces.")
    print("The final surface Sigma' is the union of Sigma and a capping surface, Scap.")
    print("The genus of the final surface, g, is the sum of the genera of the components:")
    print("g = g_sigma + g_cap")
    print("-" * 20)
    
    print("Step 4: Determine the minimal genus of the capping surface.")
    print("A key theorem in topology states that a surface with an unknotted boundary")
    print("is 'unknotted'. This means that no matter how tangled Sigma appears, it can")
    print("always be capped by a surface of genus 0 (a disk).")
    
    # Minimum possible genus for the capping surface.
    g_cap_min = 0
    print(f"Therefore, the minimum possible genus for the capping surface is g_cap = {g_cap_min}.")
    print("-" * 20)

    print("Step 5: Calculate the smallest possible genus g.")
    # Calculate the final genus g.
    g_final = g_sigma + g_cap_min

    # Output the final equation with numbers, as requested.
    print(f"The smallest genus g is found by using the minimal capping genus:")
    print(f"g = {g_sigma} + {g_cap_min}")
    print(f"g = {g_final}")

solve_surface_genus()
print("<<<10>>>", file=sys.stdout)