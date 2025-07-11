# Set the destination grid point and the constraint
max_r = 4
max_u = 8
max_k = 3 # Maximum allowed consecutive steps is 3

# DP tables initialization:
# H[r][u][k]: ways to get to (r, u) ending with k right steps
# V[r][u][k]: ways to get to (r, u) ending with k up steps
# Dimensions: (max_r+1) x (max_u+1) x (max_k+1)
H = [[[0] * (max_k + 1) for _ in range(max_u + 1)] for _ in range(max_r + 1)]
V = [[[0] * (max_k + 1) for _ in range(max_u + 1)] for _ in range(max_r + 1)]

# Iterate through all grid points (r, u) to fill the DP tables
for r in range(max_r + 1):
    for u in range(max_u + 1):
        if r == 0 and u == 0:
            continue

        # Calculate ways to reach (r, u)
        
        # Ending with a RIGHT step
        if r > 0:
            # Path ends with a single right step (e.g., ...UR)
            # This must come from (r-1, u) after an UP step.
            if r == 1 and u == 0: # Base case: The path is just 'R'
                H[r][u][1] = 1
            else:
                H[r][u][1] = V[r-1][u][1] + V[r-1][u][2] + V[r-1][u][3]
            
            # Path ends with 2 right steps (e.g., ...RR)
            # This must come from (r-1, u) after 1 right step.
            H[r][u][2] = H[r-1][u][1]

            # Path ends with 3 right steps (e.g., ...RRR)
            # This must come from (r-1, u) after 2 right steps.
            H[r][u][3] = H[r-1][u][2]

        # Ending with an UP step
        if u > 0:
            # Path ends with a single up step (e.g., ...RU)
            # This must come from (r, u-1) after a RIGHT step.
            if u == 1 and r == 0: # Base case: The path is just 'U'
                V[r][u][1] = 1
            else:
                V[r][u][1] = H[r][u-1][1] + H[r][u-1][2] + H[r][u-1][3]

            # Path ends with 2 up steps (e.g., ...UU)
            # This must come from (r, u-1) after 1 up step.
            V[r][u][2] = V[r][u-1][1]
            
            # Path ends with 3 up steps (e.g., ...UUU)
            # This must come from (r, u-1) after 2 up steps.
            V[r][u][3] = V[r][u-1][2]

# The total number of ways is the sum of all valid paths ending at the destination
h1 = H[max_r][max_u][1]
h2 = H[max_r][max_u][2]
h3 = H[max_r][max_u][3]
v1 = V[max_r][max_u][1]
v2 = V[max_r][max_u][2]
v3 = V[max_r][max_u][3]
total_ways = h1 + h2 + h3 + v1 + v2 + v3

# Print the final breakdown and the total
print(f"To find the number of unique paths from A(0,0) to B(4,8) with at most 3 consecutive moves in the same direction, we sum the possible ways to arrive at the destination B(4,8):")
print(f"Ways to arrive at (4,8) ending with 1 step Right: {h1}")
print(f"Ways to arrive at (4,8) ending with 2 steps Right: {h2}")
print(f"Ways to arrive at (4,8) ending with 3 steps Right: {h3}")
print(f"Ways to arrive at (4,8) ending with 1 step Up: {v1}")
print(f"Ways to arrive at (4,8) ending with 2 steps Up: {v2}")
print(f"Ways to arrive at (4,8) ending with 3 steps Up: {v3}")
print("\nFinal calculation:")
print(f"Total paths = {h1} + {h2} + {h3} + {v1} + {v2} + {v3} = {total_ways}")
<<<254>>>