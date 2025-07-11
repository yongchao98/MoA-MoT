import numpy as np
from itertools import combinations

class Quaternion:
    def __init__(self, a=0, b=0, c=0, d=0):
        if isinstance(a, tuple):
            self.coeffs = np.array(a, dtype=float)
        elif isinstance(a, str):
            val_map = {'1': (1,0,0,0), 'i': (0,1,0,0), 'j': (0,0,1,0), 'k': (0,0,0,1)}
            val = (0,0,0,0)
            sign = 1
            term = a
            if a.startswith('-'):
                sign = -1
                term = a[1:]
            
            if term in val_map:
                val = tuple(sign * x for x in val_map[term])
            self.coeffs = np.array(val, dtype=float)
        else:
            self.coeffs = np.array([a, b, c, d], dtype=float)

    def __repr__(self):
        return f"({self.coeffs[0]},{self.coeffs[1]}i,{self.coeffs[2]}j,{self.coeffs[3]}k)"

    def __hash__(self):
        return hash(tuple(self.coeffs))

    def __eq__(self, other):
        return np.all(self.coeffs == other.coeffs)

def get_vectors_from_string(s):
    # Remove parentheses and split into vector strings
    s = s.replace("),(", ")|(")
    vector_strings = s.strip()[1:-1].split('|')
    
    vectors = []
    for vs in vector_strings:
        parts = vs.split(',')
        vectors.append(tuple(Quaternion(p) for p in parts))
    return vectors

def get_left_mult_matrix(q):
    a, b, c, d = q.coeffs
    return np.array([
        [a, -b, -c, -d],
        [b,  a, -d,  c],
        [c,  d,  a, -b],
        [d, -c,  b,  a]
    ])

def are_linearly_dependent(v1, v2, v3):
    # We check for left linear dependence: a*v1 + b*v2 + c*v3 = 0
    # This forms a 16x12 real matrix system.
    
    # Each v is a tuple of 4 quaternions
    M_v1 = np.vstack([get_left_mult_matrix(q) for q in v1])
    M_v2 = np.vstack([get_left_mult_matrix(q) for q in v2])
    M_v3 = np.vstack([get_left_mult_matrix(q) for q in v3])
    
    A = np.hstack([M_v1, M_v2, M_v3])
    
    # The vectors are linearly dependent if the null space is non-trivial.
    # This is equivalent to the rank of the matrix being less than the number of columns.
    rank = np.linalg.matrix_rank(A)
    return rank < 12

def main():
    # The full list of vectors from the problem description
    vector_list_str = "(1,1,0,0),(1,-1,0,0),(1,0,1,0),(1,0,-1,0),(1,0,0,1),(1,0,0,-1),(0,1,1,0),(0,1,-1,0),(0,1,0,1),(0,1,0,-1),(0,0,1,1),(0,0,1,-1)," \
                      "(1,i,j,k),(1,i,-j,-k),(1,-i,j,-k),(1,-i,-j,k),(1,j,i,k),(1,j,-i,-k),(1,-j,i,-k),(1,-j,-i,k),(1,i,k,j),(1,i,-k,-j),(1,-i,k,-j),(1,-i,-k,j)," \
                      "(1,k,j,i),(1,k,-j,-i),(1,-k,i,-j),(1,-k,-i,j),(1,j,k,i),(1,j,-k,-i),(1,-j,k,-i),(1,-j,-k,i),(1,k,i,j),(1,k,-i,-j),(1,-k,i,-j),(1,-k,-i,j)"

    all_vectors_list = get_vectors_from_string(vector_list_str)
    
    # The problem has duplicate vectors in the list. We need to work with the set of unique hyperplanes.
    unique_vectors = sorted(list(set(all_vectors_list)), key=lambda v: str(v))
    
    num_unique_vectors = len(unique_vectors)
    total_triplets = len(list(combinations(unique_vectors, 3)))
    
    dependent_triplets_count = 0
    for v_triplet in combinations(unique_vectors, 3):
        if are_linearly_dependent(v_triplet[0], v_triplet[1], v_triplet[2]):
            dependent_triplets_count += 1
            
    independent_triplets_count = total_triplets - dependent_triplets_count
    
    print(f"Number of unique vectors: {num_unique_vectors}")
    print(f"Total number of triplets of vectors: C({num_unique_vectors}, 3) = {total_triplets}")
    print(f"Number of linearly dependent triplets: {dependent_triplets_count}")
    print("Dimension of H^9(M, Q) = (Total triplets) - (Dependent triplets)")
    print(f"{total_triplets} - {dependent_triplets_count} = {independent_triplets_count}")

main()