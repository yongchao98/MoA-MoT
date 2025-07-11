def solve_squarenorm_width():
    """
    This script explains the derivation of the minimum hidden-layer width
    for a shallow neural network with GeLU activations to compute the
    squared norm of an N-dimensional input vector.
    """
    print("Problem: Find the minimum hidden-layer width (H) of a shallow neural network to compute f(x) = ||x||^2 for an N-dimensional input x.")
    print("Activation function: GeLU for all hidden neurons.\n")

    print("Step 1: Decompose the target function.")
    print("The function f(x) = ||x||^2 can be written as the sum of squares of its components:")
    print("f(x) = x_1^2 + x_2^2 + ... + x_N^2")
    print("This means we can construct the network by finding a way to compute the square of a single variable, g(z) = z^2, and applying it to each input coordinate.\n")

    print("Step 2: Approximate the square function g(z) = z^2 using GeLU.")
    print("The function g(z) = z^2 is an 'even' function, since g(z) = g(-z).")
    print("The GeLU activation function is not even. To approximate an even function, we can combine GeLU(z) and GeLU(-z).")
    print("Let's analyze the function h(z) = GeLU(z) + GeLU(-z).")
    print("The Taylor expansion of GeLU(z) around z=0 is: GeLU(z) ≈ 0.5*z + 0.1995*z^2 + ...")
    print("By adding GeLU(z) and GeLU(-z), the odd-powered terms (like 0.5*z) cancel out:")
    print("h(z) = GeLU(z) + GeLU(-z) ≈ (0.5*z + 0.1995*z^2) + (-0.5*z + 0.1995*z^2) = 0.399*z^2.")
    print("Thus, z^2 can be approximated by a constant multiple of (GeLU(z) + GeLU(-z)).")
    print("This approximation for z^2 requires 2 neurons: one for GeLU(z) and one for GeLU(-z).\n")

    print("Step 3: Construct the network for f(x) = sum(x_i^2).")
    print("To compute the full sum, we apply the 2-neuron approximation for each coordinate x_i:")
    print("f(x) ≈ C * sum_{i=1 to N} [GeLU(x_i/a) + GeLU(-x_i/a)]")
    print("(Scaling by 'a' allows for arbitrary precision on a given range).\n")
    
    print("Step 4: Count the required neurons.")
    print("For each coordinate x_i, we need two neurons in the hidden layer:")
    print(" - One neuron with a weight vector that picks out x_i (e.g., w_i = (0,...,1/a,...,0)).")
    print(" - A second neuron with the negative weight vector (w_i' = (0,...,-1/a,...,0)).")
    print("For N dimensions, the weight vectors for each pair {w_i, w_i'} are unique.")
    print("This results in a total of N pairs of neurons.\n")
    
    print("Final Calculation:")
    num_neurons_per_dim = 2
    num_dimensions = "N"
    print(f"Each of the {num_dimensions} dimensions requires {num_neurons_per_dim} neurons.")
    print(f"Total minimum hidden-layer width = {num_neurons_per_dim} * {num_dimensions}")

solve_squarenorm_width()
<<<2*N>>>