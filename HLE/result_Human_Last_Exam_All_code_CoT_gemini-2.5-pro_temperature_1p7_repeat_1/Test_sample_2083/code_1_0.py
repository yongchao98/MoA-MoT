import math

def solve_shallow_network_width():
    """
    This function explains and calculates the minimum hidden-layer width
    for a shallow neural network with GeLU activation to compute the
    squared norm of an N-dimensional vector.
    """

    N_symbol = "N"

    print("### Problem Analysis ###")
    print(f"We want to find the minimum hidden-layer width (H) for a shallow neural network to compute the squared norm of an {N_symbol}-dimensional input vector x.")
    print(f"The target function is: f(x) = ||x||^2 = x_1^2 + x_2^2 + ... + x_{N_symbol}^2\n")

    print("### Step 1: Decompose the Function ###")
    print("The squared norm is a 'separable function'. It's a sum of N identical functions, each acting on a single input variable:")
    print(f"f(x) = sum(g(x_i) for i=1 to {N_symbol}), where the univariate function is g(z) = z^2.")
    print("This means we can find the minimum number of neurons 'k' to approximate g(z)=z^2, and the total width H will be N * k.\n")

    print("### Step 2: Minimum Neurons 'k' for g(z) = z^2 ###")
    print("The function g(z) = z^2 is an 'even' function, because g(z) = g(-z).")
    print("The GeLU activation function is NOT an even function. A single hidden neuron computes GeLU(w*z + b), which is also not even.")
    print("Therefore, a single neuron (k=1) cannot approximate an even function with arbitrary precision. We need k > 1.\n")

    print("### Step 3: Construction for k=2 ###")
    print("We can construct an even function by symmetrically combining two GeLU neurons:")
    print("Consider the function h(z) = GeLU(w*z) + GeLU(-w*z). This function is even.")
    print("Using a Taylor Series expansion for GeLU around z=0, we have GeLU(z) = (1/2)z + (1/sqrt(2*pi))z^2 + ...")
    print("So, h(z) = [ (1/2)wz + (1/sqrt(2*pi))(wz)^2 + ... ] + [ -(1/2)wz + (1/sqrt(2*pi))(wz)^2 - ... ]")
    print("The odd-power terms (like z) cancel out, leaving h(z) ≈ 2/sqrt(2*pi) * w^2 * z^2 for small z.")
    print("This shows a combination of two neurons can approximate z^2. Since k=1 is not enough, the minimum 'k' is 2.\n")
    k = 2

    print("### Step 4: Final Calculation ###")
    print(f"Each of the {N_symbol} input dimensions requires k={k} neurons to compute its square.")
    print("Therefore, the total minimum hidden-layer width is the product of the number of dimensions and the neurons per dimension.")
    
    # Final equation and its components
    print("\n---")
    print("Final Equation:")
    number_in_equation = k
    variable_in_equation = N_symbol
    print(f"Minimum Hidden Width = {number_in_equation} * {variable_in_equation}")
    print("---")


# Run the explanation and calculation
solve_shallow_network_width()
