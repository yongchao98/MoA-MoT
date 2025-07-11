import sys

def solve():
    """
    Calculates the maximal chromatic number of a graph G which is the
    sum of three cycles of length n.
    """

    # Step 1: Explain the interpretation of "sum" as graph join.
    print("Step 1: Interpreting the 'sum' of graphs as the graph join operation.")
    print("The graph G is taken to be the join of three cycles, G = C_n + C_n + C_n.\n")

    # Step 2: State the formula for the chromatic number of a join.
    print("Step 2: Use the formula for the chromatic number of a join:")
    print("χ(G1 + G2) = χ(G1) + χ(G2)")
    print("This implies χ(G) = χ(C_n) + χ(C_n) + χ(C_n) = 3 * χ(C_n).\n")

    # Step 3: Define a function for the chromatic number of a cycle.
    def chromatic_number_of_cycle(n):
        if n < 3:
            # Cycles are defined for n >= 3.
            return 0
        if n % 2 == 0:
            return 2
        else:
            return 3

    # Step 4: Calculate χ(G) for both even and odd n.
    print("Step 4: Calculate χ(G) based on the parity of n.")
    
    # Case 1: n is even (e.g., n=4)
    n_even = 4
    chi_cn_even = chromatic_number_of_cycle(n_even)
    chi_g_even = 3 * chi_cn_even
    print(f"For even n (e.g., n={n_even}):")
    print(f"  - χ(C_{n_even}) = {chi_cn_even}")
    print(f"  - χ(G) = 3 * {chi_cn_even} = {chi_g_even}\n")

    # Case 2: n is odd (e.g., n=5)
    n_odd = 5
    chi_cn_odd = chromatic_number_of_cycle(n_odd)
    chi_g_odd = 3 * chi_cn_odd
    print(f"For odd n (e.g., n={n_odd}):")
    print(f"  - χ(C_{n_odd}) = {chi_cn_odd}")
    print(f"  - χ(G) = 3 * {chi_cn_odd} = {chi_g_odd}\n")

    # Step 5: Determine the maximal value.
    print("Step 5: Determine the maximal chromatic number.")
    possible_values = [chi_g_even, chi_g_odd]
    maximal_chromatic_number = max(possible_values)

    # Format the final equation as requested
    numbers_in_equation = ", ".join(map(str, sorted(possible_values)))
    final_equation = f"Maximal χ(G) = max({numbers_in_equation})"
    
    print(f"The possible values for χ(G) are {set(possible_values)}.")
    print(f"The final calculation is: {final_equation} = {maximal_chromatic_number}\n")
    
    print(f"The maximal chromatic number of G is {maximal_chromatic_number}.")

solve()
<<<9>>>