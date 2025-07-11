import sympy

def solve():
    """
    This function derives the lower bound for m.
    """
    N, q, d_prime = sympy.symbols('N q d\'', positive=True, integer=True)
    m = sympy.Symbol('m', positive=True, integer=True)
    epsilon = sympy.sympify("1/(2*q)")

    print("Step 1: Define a 'hard' family of inputs.")
    print("Let v be a unit vector in R^d'. Let s be a binary vector in {-1, 1}^N.")
    print("We construct 2^N input matrices X_s, where for each matrix, the component vectors z_k are defined as z_k = s_k * v for k=1...N.")
    print("We fix the y_i components for all matrices such that for each k in {1..N}, there is an i where k is in y_i.")
    print("This construction is valid given the problem constraints (q < d', q <= N/2).\n")

    print("Step 2: Analyze the qSA function for this family of inputs.")
    print("Let s and s' be two vectors in {-1, 1}^N that differ only at index k.")
    print("s_k = 1, s'_k = -1.")
    print("Let i be an index such that k is in y_i. Such an i exists by our construction.")
    print("The qSA output for row i is:")
    print("qSA(X_s)_i = (1/q) * sum_{j in y_i} z_j = (1/q) * sum_{j in y_i} s_j * v")
    print("qSA(X_s')_i = (1/q) * sum_{j in y_i} s'_j * v\n")

    print("Step 3: Calculate the distance between the outputs.")
    print("The difference is qSA(X_s)_i - qSA(X_s')_i = (v/q) * (s_k - s'_k).")
    print("The L2 norm of this difference is ||(v/q) * (s_k - s'_k)||_2 = |s_k - s'_k|/q * ||v||_2.")
    print(f"Since s_k=1, s'_k=-1, and ||v||_2=1, the norm is |1 - (-1)|/q = 2/q.")
    dist = sympy.sympify("2/q")
    print(f"So, ||qSA(X_s)_i - qSA(X_s')_i||_2 = {dist}\n")

    print("Step 4: Use the approximation condition.")
    print("The network f must epsilon-approximate qSA, where epsilon = 1/(2*q).")
    print("By the triangle inequality, ||qSA(X_s)_i - qSA(X_s')_i|| <= ||f(X_s)_i - qSA(X_s)_i|| + ||f(X_s')_i - qSA(X_s')_i||.")
    print(f"If f(X_s) = f(X_s'), this leads to ||qSA(X_s)_i - qSA(X_s')_i|| <= 2*epsilon.")
    print(f"Substituting the values, we get {dist} <= 2 * {epsilon}, which simplifies to {dist} <= {1/q}.")
    print("This inequality (2/q <= 1/q) is false for any positive q. Thus, f(X_s) cannot equal f(X_s').\n")

    print("Step 5: Relate to the network architecture.")
    print("f(x) = g(Wx). If f(X_s) != f(X_s'), it must be that the inputs to g are different.")
    print("Therefore, Wx_s != Wx_s' for any pair s != s'.")
    print("This means the linear transformation W must map all 2^N input vectors x_s to distinct vectors in the m-dimensional hidden space.\n")

    print("Step 6: The dimensionality argument.")
    print("The input vectors x_s can be written as x_s = x_base + sum_{k=1 to N} s_k * delta_k, where delta_k are N linearly independent vectors.")
    print("These x_s vectors live in an N-dimensional affine subspace of the full input space R^(Nd).")
    print("For a linear map W to be injective on this N-dimensional space, the dimension of its image must be at least N.")
    print("The image space is R^m. Therefore, its dimension m must be at least N.\n")

    print("Step 7: Final Conclusion.")
    print(f"The analysis shows that {m} must be greater than or equal to {N}.")
    print("Asymptotically, the lower bound for m is Omega(N).")


solve()