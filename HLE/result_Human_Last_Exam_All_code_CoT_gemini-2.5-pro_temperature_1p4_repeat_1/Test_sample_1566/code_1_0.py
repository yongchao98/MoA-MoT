def solve_topology_problem():
    """
    This script solves a topology problem by logical deduction.
    It determines the number of topologically distinct continua with two specific properties
    and prints the reasoning.
    """

    # Key numbers from the problem statement that guide the logic.
    # Property 1: Number of end points > 1.
    property1_min_endpoints = 1
    # Property 2: Number of orbits = 2.
    property2_num_orbits = 2

    # The step-by-step logical deduction is presented as the output of this script.
    print("--- Step-by-Step Deduction ---")

    print("\nStep 1: Analyzing the definition of an 'end point'.")
    print("The definition states that for an end point `a`, the entire space X can be covered by a chain of open sets U_1, ..., U_N starting with `a` in U_1.")
    print("A continuum with this property is known as an 'arc-like' or 'chainable' continuum.")
    print("A fundamental theorem of topology states that any arc-like continuum has at most 2 end points.")
    print("Deduction: The number of end points, |E|, must be less than or equal to 2.")

    print("\nStep 2: Applying Property (1) to find the exact number of end points.")
    print(f"Property (1) states that |E| > {property1_min_endpoints}.")
    print(f"Combining the results of Step 1 and Property (1), we have the inequality: {property1_min_endpoints} < |E| <= {property2_num_orbits}.")
    print("Since |E| must be an integer, this forces the number of end points to be exactly 2.")
    exact_num_endpoints = 2
    print(f"Result: |E| = {exact_num_endpoints}.")

    print("\nStep 3: Identifying the unique topological structure.")
    print(f"From our deduction, the continuum X must be an arc-like continuum with exactly {exact_num_endpoints} end points.")
    print("A major theorem in topology (by R.H. Bing) states that any such continuum is homeomorphic to the closed interval [0, 1].")
    print("Therefore, any space X satisfying the conditions must be topologically equivalent to [0, 1].")

    print("\nStep 4: Verifying the structure with Property (2).")
    print(f"We must check if the interval X = [0, 1] satisfies Property (2), i.e., has exactly {property2_num_orbits} orbits under its auto-homeomorphism group.")
    print("The auto-homeomorphisms of [0, 1] are the continuous, strictly monotonic functions mapping [0, 1] to itself.")
    print(" - Orbit 1: The set of end points {0, 1}. The map h(x) = 1 - x shows that 0 and 1 are in the same orbit.")
    print(" - Orbit 2: The set of interior points (0, 1). For any two points p, q in (0, 1), a homeomorphism can be constructed to map p to q.")
    print(f"Thus, [0, 1] has exactly {property2_num_orbits} orbits, satisfying Property (2).")

    print("\nStep 5: Final Conclusion.")
    print("The analysis confirms that any continuum with the given properties must be homeomorphic to the interval [0, 1].")
    print("The question asks for the number of topologically distinct continua.")
    print("Since all such continua are topologically the same, there is only one such type.")

    final_answer = 1
    print("\n--- Final Answer ---")
    print("The properties given are:")
    print(f"(1) Number of end points is more than {property1_min_endpoints} and finite.")
    print(f"(2) Number of orbits is exactly {property2_num_orbits}.")
    print("The number of topologically distinct continua with these properties is calculated as:")
    print(f"Number of topologically distinct continua = {final_answer}")

# Run the solver.
solve_topology_problem()