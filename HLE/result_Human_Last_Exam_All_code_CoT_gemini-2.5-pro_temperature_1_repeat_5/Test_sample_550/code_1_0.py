import numpy as np
from itertools import combinations

#
# Step 1: Define quaternion arithmetic and vector operations
#
def q_mult(q1, q2):
    a1, b1, c1, d1 = q1
    a2, b2, c2, d2 = q2
    return np.array([
        a1*a2 - b1*b2 - c1*c2 - d1*d2,
        a1*b2 + b1*a2 + c1*d2 - d1*c2,
        a1*c2 - b1*d2 + c1*a2 + d1*b2,
        a1*d2 + b1*c2 - c1*b2 + d1*a2
    ], dtype=np.float64)

# Define the basis quaternions
Q_ONE = np.array([1., 0., 0., 0.])
Q_I = np.array([0., 1., 0., 0.])
Q_J = np.array([0., 0., 1., 0.])
Q_K = np.array([0., 0., 0., 1.])
Q_ZERO = np.array([0., 0., 0., 0.])

def left_mult_as_real_matrix(q):
    """Returns the 4x4 real matrix for left multiplication by q."""
    a, b, c, d = q
    return np.array([
        [a, -b, -c, -d],
        [b,  a, -d,  c],
        [c,  d,  a, -b],
        [d, -c,  b,  a]
    ])

def h_vec_to_r_vec(hv):
    """Converts a vector in H^n to a vector in R^(4n)."""
    return np.concatenate([v for v in hv])

def get_h_rank(h_vectors):
    """
    Calculates the rank of a set of vectors over the quaternions H.
    The left H-span of {v_1, ..., v_k} corresponds to the R-span of
    the block vectors { [v_s], [i*v_s], [j*v_s], [k*v_s] }_s=1..k
    where each [q*v] is the R^16 representation of the H^4 vector.
    """
    if not h_vectors:
        return 0
    
    num_h_vecs = len(h_vectors)
    n_h_dim = len(h_vectors[0])
    
    real_vectors = []
    for v_h in h_vectors: # v_h is a vector in H^4
        # v_h is a list of 4 quaternions (4-element numpy arrays)
        # To get the H-span, we form the R-span of {v, i*v, j*v, k*v}
        
        # v_r is the R^16 representation of v_h
        v_r = h_vec_to_r_vec(v_h)
        real_vectors.append(v_r)
        
        # Now find i*v_h, j*v_h, k*v_h
        iv_h = [q_mult(Q_I, q) for q in v_h]
        jv_h = [q_mult(Q_J, q) for q in v_h]
        kv_h = [q_mult(Q_K, q) for q in v_h]
        
        real_vectors.append(h_vec_to_r_vec(iv_h))
        real_vectors.append(h_vec_to_r_vec(jv_h))
        real_vectors.append(h_vec_to_r_vec(kv_h))

    # The H-rank is 1/4 of the R-rank of this expanded set
    rank_r = np.linalg.matrix_rank(np.array(real_vectors))
    return rank_r // 4

#
# Step 2: Define the set of vectors V
#
V_list = [
    (Q_ONE, Q_ONE, Q_ZERO, Q_ZERO), (Q_ONE, -Q_ONE, Q_ZERO, Q_ZERO),
    (Q_ONE, Q_ZERO, Q_ONE, Q_ZERO), (Q_ONE, Q_ZERO, -Q_ONE, Q_ZERO),
    (Q_ONE, Q_ZERO, Q_ZERO, Q_ONE), (Q_ONE, Q_ZERO, Q_ZERO, -Q_ONE),
    (Q_ZERO, Q_ONE, Q_ONE, Q_ZERO), (Q_ZERO, Q_ONE, -Q_ONE, Q_ZERO),
    (Q_ZERO, Q_ONE, Q_ZERO, Q_ONE), (Q_ZERO, Q_ONE, Q_ZERO, -Q_ONE),
    (Q_ZERO, Q_ZERO, Q_ONE, Q_ONE), (Q_ZERO, Q_ZERO, Q_ONE, -Q_ONE),
    (Q_ONE, Q_I, Q_J, Q_K), (Q_ONE, Q_I, -Q_J, -Q_K),
    (Q_ONE, -Q_I, Q_J, -Q_K), (Q_ONE, -Q_I, -Q_J, Q_K),
    (Q_ONE, Q_J, Q_I, Q_K), (Q_ONE, Q_J, -Q_I, -Q_K),
    (Q_ONE, -Q_J, Q_I, -Q_K), (Q_ONE, -Q_J, -Q_I, Q_K),
    (Q_ONE, Q_I, Q_K, Q_J), (Q_ONE, Q_I, -Q_K, -Q_J),
    (Q_ONE, -Q_I, Q_K, -Q_J), (Q_ONE, -Q_I, -Q_K, Q_J),
    (Q_ONE, Q_K, Q_J, Q_I), (Q_ONE, Q_K, -Q_J, -Q_I),
    (Q_ONE, -Q_K, Q_I, -Q_J), (Q_ONE, -Q_K, -Q_I, Q_J),
    (Q_ONE, Q_J, Q_K, Q_I), (Q_ONE, Q_J, -Q_K, -Q_I),
    (Q_ONE, -Q_J, Q_K, -Q_I), (Q_ONE, -Q_J, -Q_K, Q_I),
    (Q_ONE, Q_K, Q_I, Q_J), (Q_ONE, Q_K, -Q_I, -Q_J),
    (Q_ONE, -Q_K, Q_I, -Q_J), (Q_ONE, -Q_K, -Q_I, Q_J)
]

# Remove duplicates by converting to tuples of tuples
V_tuples = [tuple(tuple(q) for q in v) for v in V_list]
V_unique_tuples = sorted(list(set(V_tuples)))
V = [np.array(v_tuple) for v_tuple in V_unique_tuples]

#
# Step 3: Iterate through subsets to find and count circuits
#
circuit_count = 0
# We are looking for circuits of size 4 that span a 3D subspace.
subsets_of_size_4 = combinations(V, 4)

for h_vec_subset in subsets_of_size_4:
    # Check if the set of 4 vectors has rank 3
    if get_h_rank(list(h_vec_subset)) == 3:
        # Check if it's a minimal dependency (a circuit)
        # i.e., every subset of size 3 has rank 3
        is_minimal = True
        subsets_of_size_3 = combinations(h_vec_subset, 3)
        for sub_subset in subsets_of_size_3:
            if get_h_rank(list(sub_subset)) < 3:
                is_minimal = False
                break
        if is_minimal:
            circuit_count += 1

print("The dimension of the ninth cohomology group is:")
print(f"dim H^9(M, Q) = {circuit_count}")