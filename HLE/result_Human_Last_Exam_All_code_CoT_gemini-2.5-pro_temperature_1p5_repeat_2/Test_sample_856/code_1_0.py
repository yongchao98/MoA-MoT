# Step 5: Determine the number of possible structures for the odd-order part H.
# Based on the theory of metacyclic groups with cyclic Sylow subgroups, there are 2 possibilities.
num_H = 2

# Step 6: Determine the number of possible structures for the Sylow 2-subgroup P_2.
# The group must be either cyclic or generalized quaternion.
num_P2 = 2

# Step 7: Calculate the total number of manifolds.
# The total number of groups is the product of the number of choices for H and P_2.
total_manifolds = num_H * num_P2

# The problem asks to output each number in the final equation.
print(f"{num_H} * {num_P2} = {total_manifolds}")
