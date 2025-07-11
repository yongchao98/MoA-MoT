# The dimension of the system of differential equations.
n = 202000

# The number of linearly independent boundary conditions.
# Conditions given are:
# 1. x_1(T) - x_1(0) = 0
# 2. x_2(T) - x_2(0) = 0
# 3. 5(x_2(T) - x_2(0)) = 0 (redundant)
# 4. 100(x_2(T) - x_2(0)) = 0 (redundant)
# 5. -1000(x_2(T) - x_2(0)) = 0 (redundant)
# 6. 100(x_2024(T) - x_2024(0)) = 0
# Counting the unique conditions gives k=3.
k = 3

# The index of the boundary-value problem is n - k.
index = n - k

# Print the calculation and the final result.
print(f"{n} - {k} = {index}")