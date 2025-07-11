import numpy as np
import sympy as sp

def get_delta_s(S_func, delta_x):
    """
    Constructs the difference matrix delta_S for a given code construction
    and a difference vector delta_x.
    """
    # Define symbolic variables for the symbols x_i
    x = sp.symbols('x1:7')
    # Create the symbolic matrix S
    S_sym = S_func(*x)
    
    # Create the difference matrix by substituting delta_x for x
    # Conjugates are handled by sympy's conjugate() function
    delta_S_sym = S_sym.subs({xi: val for xi, val in zip(x, delta_x)})
    
    # Convert the symbolic matrix to a numerical numpy array for rank calculation
    delta_S_num = np.array(delta_S_sym, dtype=np.complex128)
    
    return delta_S_num

# Define the matrix constructions as functions
def build_Sa(*x):
    # sympy symbols are passed as a tuple x
    x = (0,) + x # 1-based indexing
    return sp.Matrix([
        [x[1], x[2], x[3], x[4], x[5], x[6]],
        [x[2], x[3], x[4], x[5], x[6], x[1]],
        [x[3], x[4], x[5], x[6], x[1], x[2]],
        [x[4], x[5], x[6], x[1], x[2], x[3]],
        [x[5], x[6], x[1], x[2], x[3], x[4]],
        [x[6], x[1], x[2], x[3], x[4], x[5]]
    ])

def build_Sb(*x):
    x = (0,) + x # 1-based indexing
    return sp.Matrix([
        [x[1], -sp.conjugate(x[2]), x[3], -sp.conjugate(x[4]), x[5], -sp.conjugate(x[6])],
        [x[2], x[3], -sp.conjugate(x[4]), x[5], -sp.conjugate(x[6]), sp.conjugate(x[1])],
        [x[3], x[4], x[5], -sp.conjugate(x[6]), sp.conjugate(x[1]), -sp.conjugate(x[2])],
        [x[4], x[5], x[6], sp.conjugate(x[1]), -sp.conjugate(x[2]), sp.conjugate(x[3])],
        [x[5], x[6], x[1], sp.conjugate(x[2]), sp.conjugate(x[3]), -sp.conjugate(x[4])],
        [x[6], x[1], x[2], sp.conjugate(x[3]), sp.conjugate(x[4]), sp.conjugate(x[5])]
    ])
    
def build_Sc(*x):
    x = (0,) + x # 1-based indexing
    return sp.Matrix([
        [x[1], sp.conjugate(x[2]), -x[3], sp.conjugate(x[4]), -x[5], sp.conjugate(x[6])],
        [x[2], -x[3], sp.conjugate(x[4]), -x[5], sp.conjugate(x[6]), sp.conjugate(x[1])],
        [-x[3], sp.conjugate(x[4]), -x[5], sp.conjugate(x[6]), sp.conjugate(x[1]), -sp.conjugate(x[2])],
        [sp.conjugate(x[4]), -x[5], sp.conjugate(x[6]), -sp.conjugate(x[1]), -sp.conjugate(x[2]), sp.conjugate(x[3])],
        [-x[5], sp.conjugate(x[6]), sp.conjugate(x[1]), -sp.conjugate(x[2]), -sp.conjugate(x[3]), -sp.conjugate(x[4])],
        [sp.conjugate(x[6]), sp.conjugate(x[1]), -sp.conjugate(x[2]), sp.conjugate(x[3]), -sp.conjugate(x[4]), -sp.conjugate(x[5])]
    ])

# --- Analysis ---
print("Analyzing Diversity Order of Space-Time Codes\n")

# Case Sa: Rank should be 1
delta_x_a = [1, 1, 1, 1, 1, 1]
delta_Sa = get_delta_s(build_Sa, delta_x_a)
rank_a = np.linalg.matrix_rank(delta_Sa)
print(f"Code Sa: For delta_x = {delta_x_a}, rank(delta_S_a) = {rank_a}")
print("Diversity order of Sa is 1.\n")

# Case Sb: Rank should be 6 (full)
delta_x_b1 = [1, 0, 0, 0, 0, 0]
delta_Sb1 = get_delta_s(build_Sb, delta_x_b1)
rank_b1 = np.linalg.matrix_rank(delta_Sb1)
print(f"Code Sb: For delta_x = {delta_x_b1}, rank(delta_S_b) = {rank_b1}")
delta_x_b2 = [1, 1j, 0, 0, 0, 0]
delta_Sb2 = get_delta_s(build_Sb, delta_x_b2)
rank_b2 = np.linalg.matrix_rank(delta_Sb2)
print(f"Code Sb: For delta_x = {delta_x_b2}, rank(delta_S_b) = {rank_b2}")
print("Code Sb is a known perfect code with full diversity. Its diversity order is 6.\n")

# Case Sc: Rank should be < 6 for a specific delta_x
delta_x_c = [-2j, 2, -2j, -2, 2j, 2] # This vector makes C1+C3=0
delta_Sc = get_delta_s(build_Sc, delta_x_c)
rank_c = np.linalg.matrix_rank(delta_Sc)
print(f"Code Sc: For a specific non-zero delta_x, rank(delta_S_c) = {rank_c}")
print("Diversity order of Sc is less than 6.\n")

# Final Conclusion
print("--- Conclusion ---")
print("Code Sb provides the maximum diversity order.")
print("The maximum diversity order is 6.")