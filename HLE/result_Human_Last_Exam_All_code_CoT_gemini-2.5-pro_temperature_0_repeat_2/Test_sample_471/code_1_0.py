def solve_torus_critical_points():
    """
    Calculates the minimal number of critical points for a smooth function
    on a 2-torus (T^2) using principles from Morse theory.
    """
    print("To find the minimal number of critical points for a smooth function f: T^2 -> R, we use Morse theory.")
    print("This theory relates the number of critical points to the topology of the domain, the 2-torus.")
    print("\nStep 1: Define the Betti numbers of the 2-torus (T^2).")
    print("Betti numbers describe the topological 'holes' of a space.")

    # Betti numbers for the 2-torus
    b0 = 1  # Number of connected components
    b1 = 2  # Number of one-dimensional 'circular' holes (the two fundamental loops of a torus)
    b2 = 1  # Number of two-dimensional 'voids' (the surface itself)

    print(f" - b_0 = {b0} (The torus is one connected piece)")
    print(f" - b_1 = {b1} (The torus has two independent loops, like latitude and longitude)")
    print(f" - b_2 = {b2} (The torus encloses one surface)")

    print("\nStep 2: Apply the Morse Inequalities.")
    print("The inequalities state that the number of critical points of index k (c_k) must be at least the k-th Betti number (b_k).")
    print(" - c_0 (minima) >= b_0")
    print(" - c_1 (saddles) >= b_1")
    print(" - c_2 (maxima) >= b_2")

    # The minimal number of critical points of each type is given by the Betti numbers.
    min_c0 = b0
    min_c1 = b1
    min_c2 = b2

    print("\nStep 3: Calculate the minimal number of each type of critical point.")
    print(f" - Minimal number of minima (c_0) = {min_c0}")
    print(f" - Minimal number of saddle points (c_1) = {min_c1}")
    print(f" - Minimal number of maxima (c_2) = {min_c2}")

    # The total minimal number of critical points is the sum of these minima.
    # This lower bound is achievable, for example, by the height function on a standard
    # torus standing upright in R^3.
    total_min_critical_points = min_c0 + min_c1 + min_c2

    print("\nStep 4: Sum the minima to find the total minimal number of critical points.")
    print("The final equation is the sum of the minimal number of minima, saddles, and maxima.")
    print("\n--- Final Calculation ---")
    print(f"{min_c0} + {min_c1} + {min_c2} = {total_min_critical_points}")
    print("-------------------------")

    print(f"\nConclusion: Any smooth function on the 2-torus must have at least {total_min_critical_points} critical points.")

# Execute the function to print the solution
solve_torus_critical_points()