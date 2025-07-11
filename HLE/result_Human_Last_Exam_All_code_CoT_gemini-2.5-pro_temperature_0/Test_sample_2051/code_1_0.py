def solve_moduli_volume_problem():
    """
    Solves the two-part question about the volume of moduli spaces of oriented ribbon graphs.

    Part (a) is a theoretical question about continuity.
    Part (b) requires calculating the degree of a specific polynomial volume.
    """

    # Part (a): Continuity
    # The volume of a geometric object like a moduli space is a continuous function.
    # While being "piecewise polynomial" does not in general guarantee continuity,
    # for the specific case of these Weil-Petersson volumes, the function is continuous
    # across the boundaries of the polynomial pieces (cells).
    answer_a = "Yes"

    # Part (b): Degree of the polynomial Z_{0,3,1}
    # The degree of homogeneity of the volume polynomial Z_{g,n} is given by the formula:
    # Degree = 2*g - 2 + n
    # where n = n_+ + n_-.

    # Given values for the specific case
    g = 0
    n_plus = 3
    n_minus = 1

    # Calculate the total number of boundaries, n
    n = n_plus + n_minus

    # Calculate the degree of the polynomial
    degree = 2 * g - 2 + n
    answer_b = degree

    # --- Output ---
    print("This script calculates the degree of the polynomial Z_{g, n_+, n_-}(L_+ | L_-) for the specified case.")
    print("\n--- Calculation for Part (b) ---")
    print("The formula for the degree is: Degree = 2*g - 2 + n")
    print(f"Given values are: g = {g}, n_+ = {n_plus}, n_- = {n_minus}")
    print(f"First, we find the total number of boundaries, n:")
    print(f"n = n_+ + n_- = {n_plus} + {n_minus} = {n}")
    print("\nNow, we substitute the values into the degree formula:")
    # The final code needs to output each number in the final equation
    print(f"Degree = 2 * {g} - 2 + {n}")
    print(f"Degree = {2 * g} - 2 + {n}")
    print(f"Degree = {degree}")

    print("\n--- Final Answer ---")
    # Print the final answer in the required format
    print(f"(a) {answer_a}; (b) {answer_b}")


solve_moduli_volume_problem()