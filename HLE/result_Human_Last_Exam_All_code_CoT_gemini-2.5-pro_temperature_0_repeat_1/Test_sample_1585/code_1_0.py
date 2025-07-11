# Step 1: Define the hypothesized values for n_2 and n_3 based on the reasoning.
# n_2 is the smallest number of vertices for the 2-planar graph.
# Based on the structural hypothesis n_k = m * k with a base size m=10, n_2 = 10 * 2.
n_2 = 20

# n_3 is the smallest number of vertices for the 3-planar graph.
# Based on the same structural hypothesis, n_3 = 10 * 3.
n_3 = 30

# Step 2: Calculate the required expression.
result = (n_2 + n_3) * n_2

# Step 3: Print the final equation with all the numbers.
print(f"({n_2} + {n_3}) * {n_2} = {result}")

# Step 4: Output the final answer in the specified format.
print(f"<<<{result}>>>")