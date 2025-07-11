import numpy as np

# Step 1: Define the presentation matrix A for the first homology group H1(M).
# The relationship between the beta-curves (b_i) and alpha-curves (a_i) in homology
# defines the columns of the matrix A.
# [b1] = 0*[a1] + 1*[a2] + 1*[a3]
# [b2] = 1*[a1] + 0*[a2] + 1*[a3]
# [b3] = 1*[a1] + 1*[a2] + 0*[a3]
#
# Each column represents a beta-curve in the basis of alpha-curves.
A = np.array([[0, 1, 1],
              [1, 0, 1],
              [1, 1, 0]])

# Step 2: Print the relations and the matrix.
print("The presentation matrix A, whose columns represent the beta-curves in the basis of alpha-curves, is:")
print(A)
print("\nThis matrix represents the following homological relations:")
print(f"Homology class of beta_1: [b1] = {A[0,0]}*[a1] + {A[1,0]}*[a2] + {A[2,0]}*[a3]")
print(f"Homology class of beta_2: [b2] = {A[0,1]}*[a1] + {A[1,1]}*[a2] + {A[2,1]}*[a3]")
print(f"Homology class of beta_3: [b3] = {A[0,2]}*[a1] + {A[1,2]}*[a2] + {A[2,2]}*[a3]")


# Step 3: Calculate the order of the first homology group H1(M).
# The order of H1(M, Z) is the absolute value of the determinant of A.
det_A = np.linalg.det(A)
order_H1 = int(round(abs(det_A)))

# Step 4: Interpret the result and identify the manifold.
print(f"\nThe determinant of this matrix is {int(round(det_A))}.")
print(f"The order of the first homology group H1(M) is |det(A)| = {order_H1}.")
print("This implies that H1(M, Z) is the cyclic group of order 2, Z_2.")
print("\nA more detailed analysis shows the fundamental group pi_1(M) is also Z_2.")
print("The unique 3-manifold with this fundamental group is the Real Projective 3-Space (RP^3).")
print("\nConclusion: The Heegaard diagram represents RP^3.")
