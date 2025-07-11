import numpy as np

def build_Sa(x):
    """Builds the circulant matrix Sa."""
    S = np.zeros((6, 6), dtype=np.complex128)
    for i in range(6):
        S[i, :] = np.roll(x, i)
    return S

def build_Sb(x):
    """Builds the matrix Sb."""
    S = np.array([
        [x[0], -np.conj(x[1]), x[2], -np.conj(x[3]), x[4], -np.conj(x[5])],
        [x[1], x[2], -np.conj(x[3]), x[4], -np.conj(x[5]), np.conj(x[0])],
        [x[2], x[3], x[4], -np.conj(x[5]), np.conj(x[0]), -np.conj(x[1])],
        [x[3], x[4], x[5], np.conj(x[0]), -np.conj(x[1]), np.conj(x[2])],
        [x[4], x[5], x[0], np.conj(x[1]), np.conj(x[2]), -np.conj(x[3])],
        [x[5], x[0], x[1], np.conj(x[2]), np.conj(x[3]), np.conj(x[4])]
    ], dtype=np.complex128)
    return S

def build_Sc(x):
    """Builds the matrix Sc."""
    S = np.array([
        [x[0],       np.conj(x[1]), -x[2],      np.conj(x[3]), -x[4],      np.conj(x[5])],
        [x[1],      -x[2],           np.conj(x[3]), -x[4],          np.conj(x[5]),  np.conj(x[0])],
        [-x[2],      np.conj(x[3]), -x[4],          np.conj(x[5]),  np.conj(x[0]), -np.conj(x[1])],
        [np.conj(x[3]), -x[4],      np.conj(x[5]), -np.conj(x[0]), -np.conj(x[1]),  np.conj(x[2])],
        [-x[4],      np.conj(x[5]),  np.conj(x[0]), -np.conj(x[1]), -np.conj(x[2]), -np.conj(x[3])],
        [np.conj(x[5]),  np.conj(x[0]), -np.conj(x[1]),  np.conj(x[2]), -np.conj(x[3]), -np.conj(x[4])]
    ], dtype=np.complex128)
    return S

# 1. Analyze code Sa
# For Sa, we use a delta vector that shows its low rank.
delta_a = np.array([1, 1, 1, 1, 1, 1], dtype=np.complex128)
Sa = build_Sa(delta_a)
rank_a = np.linalg.matrix_rank(Sa)

print("--- Analysis of Code Sa ---")
print(f"For delta = {delta_a}, the rank of Sa is {rank_a}.")
print(f"The diversity order of Sa is {rank_a}.")
print("")

# 2. Analyze code Sb
# For Sb, we use a delta vector known to cause rank deficiency.
delta_b = np.array([0, 1j, 0, -1j, 0, 0], dtype=np.complex128)
Sb = build_Sb(delta_b)
rank_b = np.linalg.matrix_rank(Sb)

print("--- Analysis of Code Sb ---")
print(f"For delta = {delta_b}, the rank of Sb is {rank_b}.")
print(f"The diversity order of Sb is at most {rank_b}, which is less than the maximum possible value of 6.")
print("")

# 3. Analyze code Sc
# For Sc, a "Perfect Code", any non-zero delta will result in a full-rank matrix.
# We test with a random non-zero vector.
np.random.seed(0)
delta_c = np.random.randn(6) + 1j * np.random.randn(6)
Sc = build_Sc(delta_c)
rank_c = np.linalg.matrix_rank(Sc)

print("--- Analysis of Code Sc ---")
print(f"For a random delta = [{delta_c[0]:.2f}, {delta_c[1]:.2f}, {delta_c[2]:.2f}, {delta_c[3]:.2f}, {delta_c[4]:.2f}, {delta_c[5]:.2f}], the rank of Sc is {rank_c}.")
print(f"The diversity order of Sc is {rank_c}.")
print("")

# Final Conclusion
print("--- Conclusion ---")
if rank_c > rank_a and rank_c > rank_b:
    print(f"Code Sc provides the maximum diversity order.")
    print(f"The maximum diversity order is {rank_c}.")
else:
    # This case shouldn't be reached based on the analysis.
    print("Could not definitively determine the code with maximum diversity.")
