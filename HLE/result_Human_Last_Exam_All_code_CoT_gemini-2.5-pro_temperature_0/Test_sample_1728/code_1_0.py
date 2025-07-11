def solve_maximal_chromatic_number():
    """
    This script explains and calculates the maximal chromatic number of a graph G
    that is the sum of three cycles of length n.
    """
    print("Step 1: Defining the graph G and the 'sum' operation.")
    print("The graph G is the sum of three cycles of length n, which we can write as G = C_n + C_n + C_n.")
    print("The 'sum' operation in this context is the graph join. The join of several graphs is formed by connecting every vertex of each graph to every vertex of the other graphs.\n")

    print("Step 2: Applying the chromatic number formula for graph joins.")
    print("A key theorem in graph theory states that the chromatic number of a join of graphs is the sum of their individual chromatic numbers.")
    print("Therefore, χ(G) = χ(C_n) + χ(C_n) + χ(C_n).\n")

    print("Step 3: Finding the chromatic number of a cycle, χ(C_n).")
    print("The chromatic number of a cycle C_n depends on its length n (for n ≥ 3):")
    print(" - If n is even, χ(C_n) = 2.")
    print(" - If n is odd,  χ(C_n) = 3.\n")

    print("Step 4: Maximizing the chromatic number of G.")
    print("To find the maximal chromatic number of G, we must choose n such that χ(C_n) is at its maximum.")
    max_chi_cn = 3
    print(f"The maximum value for χ(C_n) is {max_chi_cn}, which occurs when n is any odd number.\n")

    print("Step 5: Calculating the final result.")
    print("Substituting the maximal value of χ(C_n) into our equation for χ(G), we get the final equation:")
    
    # Define the values for the final equation
    chi_cycle_1 = max_chi_cn
    chi_cycle_2 = max_chi_cn
    chi_cycle_3 = max_chi_cn
    
    # Calculate the total maximal chromatic number
    max_chi_g = chi_cycle_1 + chi_cycle_2 + chi_cycle_3
    
    # Print the final equation with each number
    print(f"{chi_cycle_1} + {chi_cycle_2} + {chi_cycle_3} = {max_chi_g}")

if __name__ == "__main__":
    solve_maximal_chromatic_number()