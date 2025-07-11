import math

def solve_min_width():
    """
    This function explains the derivation for the minimum hidden-layer width
    required for a shallow neural network to compute the squared norm of an
    N-dimensional vector.
    """

    print("### Step-by-step Derivation ###")
    print("\n--- 1. The Neural Network Model ---")
    print("Let the input be an N-dimensional vector x.")
    print("The target function is the squared norm: f(x) = ||x||^2 = x_1^2 + ... + x_N^2.")
    print("The network is a shallow network with one hidden layer of width H and a linear output neuron.")
    print("The hidden neurons use the GELU activation function.")
    print("The network's output y(x) is given by:")
    print("y(x) = sum_{i=1 to H} [c_i * GELU(w_i^T * x + b_i)] + b_out")
    print("where w_i are N-dimensional weight vectors, and c_i, b_i, b_out are scalar biases and output weights.")

    print("\n--- 2. Local Approximation via Taylor Expansion ---")
    print("To approximate a quadratic function, we analyze the network's behavior for small x around x=0.")
    print("We can use a Taylor series expansion of y(x) around x=0.")
    print("For simplicity, let's set all hidden biases b_i = 0. This is a valid choice that simplifies the analysis without loss of generality for the result.")
    print("The Taylor expansion of GELU(z) around z=0 is: GELU(z) ≈ z/2 + z^2/sqrt(2*pi) + O(z^4).")
    print("Substituting z = w_i^T * x into the expansion of y(x):")
    print("y(x) ≈ sum_{i=1 to H} [c_i * ( (w_i^T*x)/2 + (w_i^T*x)^2/sqrt(2*pi) )] + b_out")
    print("y(x) ≈ b_out + (1/2 * sum(c_i*w_i))^T * x + (1/sqrt(2*pi)) * x^T * (sum(c_i*w_i*w_i^T)) * x")

    print("\n--- 3. Deriving Constraints ---")
    print("We want y(x) to approximate ||x||^2 = x^T * I * x, where I is the N x N identity matrix.")
    print("By comparing the terms of our y(x) expansion with the target function, we get the following constraints:")
    print(" a) Constant term: b_out = 0")
    print(" b) Linear term: The linear term in ||x||^2 is zero. So, we must have: sum_{i=1 to H} (c_i * w_i) = 0")
    print(" c) Quadratic term: The quadratic parts must match. This gives: sum_{i=1 to H} (c_i * w_i * w_i^T) = k * I")
    print("    where w_i*w_i^T is an outer product (an NxN matrix) and k = sqrt(2*pi) is a positive constant.")

    print("\n--- 4. Minimum Width (H) - The Lower Bound ---")
    print("Let V be the N x H matrix whose columns are the weight vectors w_1, ..., w_H.")
    print("Let c be the H x 1 column vector of the output weights c_1, ..., c_H.")
    print("The constraints in matrix form are:")
    print(" (L) V * c = 0")
    print(" (Q) V * C * V^T = k * I, where C is the HxH diagonal matrix of c_i.")
    print("\nFrom constraint (Q), the matrix V*C*V^T must have rank N (since k*I is rank N).")
    print("The rank of a product of matrices is at most the minimum of their ranks. rank(V) must be at least N.")
    print("Since V is an N x H matrix, its rank is at most min(N, H). Thus, we must have H >= N.")
    print("\nFrom constraint (L), the vector c must be in the null space of V.")
    print("For a non-trivial solution (where not all c_i are zero), the null space of V must have a dimension greater than 0.")
    print("The dimension of the null space of V is given by the rank-nullity theorem: dim(null(V)) = H - rank(V).")
    print("Since rank(V) = N, we have dim(null(V)) = H - N.")
    print("For a non-trivial c to exist, we need H - N > 0, which implies H > N.")
    print("Since H must be an integer, the minimum possible value for H is N + 1.")

    print("\n--- 5. Sufficiency of H = N + 1 ---")
    print("We need to show that H = N+1 is not just necessary, but also sufficient.")
    print("This requires finding N+1 vectors w_i in R^N and coefficients c_i that satisfy the constraints.")
    print("A known construction exists: Choose c_i = 1 for all i=1,...,N+1.")
    print("The constraints become:")
    print(" (L) sum(w_i) = 0")
    print(" (Q) sum(w_i * w_i^T) = k * I")
    print("A set of vectors satisfying these properties is known as a 'simplex frame'.")
    print("The N+1 vertices of a regular N-simplex centered at the origin form such a set of vectors.")
    print("For any N, such a set of N+1 vectors can be constructed.")
    print("Therefore, a solution exists for H = N + 1.")

    print("\n--- 6. Conclusion ---")
    print("The minimum required width H is greater than N, and a width of N+1 is sufficient.")
    print("Therefore, the minimum hidden-layer width is N+1.")
    print("\nThe final equation is H = N + 1.")
    print("To explicitly show the numbers in the final equation: H = 1 * N + 1")

if __name__ == '__main__':
    solve_min_width()
