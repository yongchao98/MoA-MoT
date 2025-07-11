# Define the sequence of characters for X1 through X10.
# We'll interpret them as the digits '1' through '9' and '0'.
X = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '0']

# Define the 1-based indices for Y and Z from the problem statement.
y_indices = [1, 2, 3, 4, 2, 5]
z_indices = [3, 6, 1, 7, 8, 9, 10, 2, 6]

# Build the lists of characters for Y and Z using the indices.
# We subtract 1 from each index to match Python's 0-based list indexing.
y_components = [X[i - 1] for i in y_indices]
z_components = [X[i - 1] for i in z_indices]

# Join the component characters to form the final strings for Y and Z.
Y = "".join(y_components)
Z = "".join(z_components)

# Fulfill the requirement to "output each number in the final equation".
print("Y =", ' '.join(y_components))
print("Z =", ' '.join(z_components))

# Print the final answer "Y Z".
print(f"\n{Y} {Z}")