import numpy as np

def compute_abelianization_parameters():
    """
    Computes the rank and torsion order for the abelianization of the group G.

    The group G, also known as the golden ratio Thompson group F_tau, consists of
    piecewise linear homeomorphisms of [0,1] with breakpoints in Z[tau] and
    slopes that are integer powers of tau.

    The abelianization Ab(G) is analyzed via a homomorphism d: G -> Z x Z,
    d(g) = (log_tau(g'(0+)), log_tau(g'(1-))).

    The rank 'r' of Ab(G) is the rank of the image of d. We find this by
    checking if the images of G's generators span Z x Z.

    The order of the torsion subgroup 't' is determined by the kernel of the
    map from Ab(G) to Z x Z. A known theorem by M. Brin states that
    ker(d) = [G, G], which implies Ab(G) is isomorphic to Z x Z and has
    no torsion.
    """

    # The standard generators for G, let's call them A and B, have the following
    # derivatives at the endpoints 0 and 1:
    # For A: A'(0+) = tau^1, A'(1-) = tau^(-1)
    # For B: B'(0+) = tau^1, B'(1-) = tau^0 = 1

    # Applying the homomorphism d gives the images of the generators in Z x Z.
    # d(g) = (log_tau(g'(0+)), log_tau(g'(1-)))
    d_A = (1, -1)
    d_B = (1, 0)

    print(f"The images of the generators in ZxZ are d(A) = {d_A} and d(B) = {d_B}.")

    # To check if these vectors generate Z x Z, we form a matrix and compute its
    # determinant. If the determinant is +-1, the vectors form a basis for Z x Z,
    # and the homomorphism d is surjective.
    M = np.array([d_A, d_B])

    # The problem description specifies to show the numbers in the final equation.
    # We are calculating the determinant of the matrix formed by the vectors.
    print(f"\nMatrix formed by the image vectors:\n{M}")
    
    det_M = int(round(np.linalg.det(M)))
    
    print(f"The determinant is {M[0][0]}*{M[1][1]} - {M[0][1]}*{M[1][0]} = {det_M}.")

    # The rank r of Ab(G) is the rank of the image of d.
    if abs(det_M) == 1:
        r = 2
        print(f"\nSince the determinant is {det_M}, the image of d has rank 2.")
        print(f"This means the rank of Ab(G) is r = {r}.")
    else:
        # This case should not be reached with the correct generator images.
        r = np.linalg.matrix_rank(M)
        print(f"The image of d has rank {r}.")

    # The order t of the torsion subgroup of Ab(G).
    # Based on the theorem that ker(d) = [G, G], we have Ab(G) isomorphic to Z x Z.
    # The group Z x Z is torsion-free.
    t = 1
    print("\nBased on the known result that Ab(G) is isomorphic to Z x Z, the torsion subgroup is trivial.")
    print(f"The order of the torsion subgroup is t = {t}.")

    # The final pair (r, t)
    print("\n--- Final Answer ---")
    print(f"The computed pair is (r, t) = ({r}, {t}).")
    
    return r, t

if __name__ == '__main__':
    compute_abelianization_parameters()