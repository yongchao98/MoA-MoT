import numpy as np

def poly_to_str(p):
    """
    Custom function to convert a numpy.poly1d object to a clean string representation
    for display. Example: '1 s + 1' becomes 's + 1'.
    """
    # Use numpy's built-in string conversion with variable 's'
    p_str = np.poly2str(p, 's')
    
    # Clean up the string for better readability
    p_str = p_str.replace('1 s', 's')
    if p_str.startswith('s'):
        p_str = ' ' + p_str
    p_str = p_str.replace('+ -', '- ')
    p_str = p_str.replace('-  ', '- ')
    
    return p_str

def print_poly_matrix(M, title=""):
    """Prints a matrix of polynomials with a title."""
    if title:
        print(title)
    
    # Find the maximum width needed for any element string
    max_len = 0
    str_matrix = np.empty(M.shape, dtype=object)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            s = poly_to_str(M[i,j])
            str_matrix[i,j] = s
            if len(s) > max_len:
                max_len = len(s)

    # Print the matrix with aligned columns
    for i in range(M.shape[0]):
        row_str = "  ".join([f"{s:^{max_len}}" for s in str_matrix[i,:]])
        print(f"[ {row_str} ]")
    print("-" * (len(row_str) + 4))


# Define the polynomial variable 's'
s = np.poly1d([1, 0])

# Define the polynomial matrices P1 and P2
P1 = np.array([
    [s**2 + s, -s],
    [-s**2 - 1, s**2]
], dtype=object)

P2 = np.array([
    [s, np.poly1d(0)],
    [-s - 1, np.poly1d(1)]
], dtype=object)

# Step 1: Stack P1 and P2 to form the initial matrix M
M = np.vstack((P1, P2))
print_poly_matrix(M, "Initial stacked matrix M = [P1; P2]:")

# Step 2: Apply elementary row operations
# The row M[3,:] = [-s - 1, 1] contains a constant '1', which is an excellent pivot.
# Let's swap it with the first row to make it the primary pivot row.
M[[0, 3], :] = M[[3, 0], :]
print_poly_matrix(M, "After swapping R1 and R4:")

# Now, use the new R1 (M[0,:]) to eliminate elements in the second column of other rows.
# R2 -> R2 - s^2 * R1
multiplier = M[1, 1]
M[1, :] = M[1, :] - multiplier * M[0, :]

# R4 -> R4 + s * R1
multiplier = M[3, 1]
M[3, :] = M[3, :] - multiplier * M[0, :]
print_poly_matrix(M, "After zeroing out column 2 using R1:")

# The matrix now has two rows with a single non-zero polynomial in the first column.
# M[1,:] = [s^3 - 1, 0] and M[2,:] = [s, 0].
# We can reduce these further using polynomial division (Euclidean algorithm).
# R2 -> R2 - q * R3, where q is the quotient of (s^3-1) / s.
p1 = M[1, 0]
p2 = M[2, 0]
quotient, remainder = np.polydiv(p1, p2)
M[1, :] = M[1, :] - quotient * M[2, :]
print_poly_matrix(M, "After reducing R2 with R3:")

# Now use the new R2 (M[1,:] = [-1, 0]) to eliminate R3 (M[2,:] = [s, 0]).
# R3 -> R3 - q * R2, where q is the quotient of s / -1.
p1 = M[2, 0]
p2 = M[1, 0]
quotient, remainder = np.polydiv(p1, p2)
M[2, :] = M[2, :] - quotient * M[1, :]
print_poly_matrix(M, "After reducing R3 with the new R2:")

# The matrix is now in row-echelon form. The top two non-zero rows form the GCRD.
G = M[:2, :]

# Step 3: Simplify the GCRD for a canonical representation.
# This is done by left-multiplication with unimodular matrices.
# 1. Multiply the second row by -1 to make the leading coefficient positive.
G[1, :] = -G[1, :]
# 2. Swap the two rows to have the identity matrix element at the top-left.
G[[0, 1], :] = G[[1, 0], :]

# Final result
print("\nThe greatest common right divisor is:")
max_len = max(len(poly_to_str(p)) for p in G.flatten())
for i in range(G.shape[0]):
    row_str = " , ".join([f"{poly_to_str(p):>{max_len}}" for p in G[i,:]])
    print(f"[{row_str}]")
