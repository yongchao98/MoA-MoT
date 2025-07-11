import math

def combinations(n, k):
    """Calculates the binomial coefficient 'n choose k'."""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def calculate_f_vector_of_neighborly_polytope(d, v):
    """
    Calculates the f-vector of a neighborly d-polytope with v vertices.
    This corresponds to the polytope that maximizes the number of faces of all dimensions.
    """
    # 1. Calculate the g-vector
    # For a neighborly d-polytope, the g-vector components are determined by g_1.
    g = [0] * (d // 2 + 1)
    g[0] = 1
    if d >= 2:
        g1 = v - (d + 1)
        g[1] = g1
    
    # For a neighborly polytope, g_k is C(g_1 + k - 1, k)
    for k in range(2, len(g)):
        g[k] = combinations(g[1] + k - 1, k)
    
    # 2. Calculate the h-vector from the g-vector
    # h_k = sum of g_i from i=0 to k
    h = [0] * (d + 1)
    current_g_sum = 0
    for k in range(len(g)):
        current_g_sum += g[k]
        h[k] = current_g_sum
        
    # Apply Dehn-Sommerville equations (h_k = h_{d-k})
    for k in range(d // 2 + 1, d + 1):
        h[k] = h[d - k]

    # 3. Calculate the f-vector from the h-vector
    f = [0] * d
    f[0] = v # Number of vertices is given
    # f_{k-1} = sum_{i=0 to k} C(d-i, k-i) * h_i
    for k in range(1, d): # Iterate for f_1, f_2, ..., f_{d-1}
        f_k_minus_1 = 0
        for i in range(k + 1):
            f_k_minus_1 += combinations(d - i, k - i) * h[i]
        f[k] = f_k_minus_1
    
    # As the f_k calculation is complex, a simpler set of relations can be used:
    f[0] = v
    f_val = [0] * d
    f_val[0] = v
    for k in range(1, d): # f_1 to f_{d-1}
        term_sum = 0
        for i in range(k + 1):
             term_sum += combinations(d - i, k - i) * h[i]
        f_val[k] = term_sum

    # A simpler set of relations to calculate f from h:
    f_calc = {}
    f_calc[-1] = 1
    for j in range(d): # f_0 to f_{d-1}
        f_calc[j] = sum(combinations(d - i, j - i + 1) * f_calc[i-1] for i in range(j + 1))
        
    # The standard relations from Ziegler "Lectures on Polytopes"
    f_vector = [0] * d
    for j in range(d): # Computes f_j
        f_vector[j] = sum(combinations(d - i, d - j) * h[i] for i in range(j + 1))
    
    # Using the sum relation from Wikipedia (h -> f)
    # f_{k-1} = sum_{i=0 to k} C(d-i, d-k) * h_i
    f_final = [-1]*d
    f_final[0] = v
    for k in range(2, d+1):
      val = sum([combinations(d-i,d-k)*h[i] for i in range(k+1)])
      f_final[k-1] = val

    return tuple(f_final)

# Parameters from the problem
d = 4  # dimension
f0 = 6 # number of vertices

# We solve for the polytope with the maximal number of faces, which must be simplicial.
f = calculate_f_vector_of_neighborly_polytope(d, f0)

print("Based on the analysis that any 4-polytope with 6 vertices must be simplicial,")
print("we calculate the f-vector for the one with the maximal number of 2-faces.")
print("This is the neighborly 4-polytope with 6 vertices.\n")

f0, f1, f2, f3 = f

print(f"The calculated f-vector is (f_0, f_1, f_2, f_3): {f}\n")
print("The components are:")
print(f"f_0 (vertices) = {f0}")
print(f"f_1 (edges) = {f1}")
print(f"f_2 (2-faces) = {f2}")
print(f"f_3 (facets/cells) = {f3}\n")

print("Checking with Euler's Polytope Formula (f_0 - f_1 + f_2 - f_3 = 0):")
result = f0 - f1 + f2 - f3
print(f"{f0} - {f1} + {f2} - {f3} = {result}")
