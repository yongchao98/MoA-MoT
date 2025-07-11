import numpy as np

# Set up the system of linear equations Ac = b based on the data from the 5 graphs.
# Each row in A corresponds to a graph, and columns correspond to the graph properties.
# The vector b contains the number of walks for each graph.

# A = [[e, k, p, S2, S3] for each graph]
A = np.array([
    [3, 1, 0, 3, 0],  # K3
    [4, 0, 4, 4, 0],  # C4
    [3, 0, 1, 2, 0],  # P4
    [3, 0, 0, 3, 1],  # K1,3
    [1, 0, 0, 0, 0]   # P2
])

# b = [Number of closed tree-like walks of length 6]
b = np.array([
    42,  # K3
    80,  # C4
    36,  # P4
    54,  # K1,3
    2    # P2
])

# Solve the system of equations for the coefficients c
try:
    c = np.linalg.solve(A, b)
    
    # The result might be in float, so we round and convert to integers.
    c_int = np.round(c).astype(int)
    
    c1, c2, c3, c4, c5 = c_int
    
    print(f"The calculated coefficients are:")
    print(f"c_1 = {c1}")
    print(f"c_2 = {c2}")
    print(f"c_3 = {c3}")
    print(f"c_4 = {c4}")
    print(f"c_5 = {c5}")
    
    print("\nThe expression for the number of closed tree-like walks of length 6 is:")
    print(f"N = {c1} * e + {c2} * k + {c3} * p + {c4} * sum(deg(v) choose 2) + {c5} * sum(deg(v) choose 3)")
    
    print("\nThe final coefficients in order c_1, c_2, c_3, c_4, c_5 are:")
    # Print the final answer in the requested format
    final_answer = ",".join(map(str, c_int))
    print(final_answer)

except np.linalg.LinAlgError:
    print("The system of equations could not be solved. The matrix might be singular.")

<<<2,0,6,12,12>>>