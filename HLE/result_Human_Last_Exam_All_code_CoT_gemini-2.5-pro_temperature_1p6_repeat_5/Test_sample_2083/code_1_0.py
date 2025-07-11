import sys

def solve_neural_network_width():
    """
    Determines and explains the minimum hidden-layer width for a shallow neural network
    to compute the squared norm of an N-dimensional vector using GeLU activation.
    """

    explanation = """
Step 1: Establishing a lower bound for the hidden width (H)

Let the network's output function be f(x), which approximates the target function g(x) = ||x||^2 = \sum_{i=1}^{N} x_i^2.
The network has one hidden layer with H neurons, so its output (before the final output activation) is of the form:
f(x) = sum_{j=1 to H} c_j * GeLU(w_j^T * x + b_j) + d

Let's consider the Hessian matrix of second derivatives.
The Hessian of the target function g(x) is H_g.
The partial derivative dg/dx_k = 2*x_k.
The second partial derivative d^2g/dx_k^2 = 2.
The mixed partial derivative d^2g/dx_k*dx_l = 0 for k != l.
So, H_g is a diagonal matrix with all 2s on the diagonal: H_g = 2*I, where I is the N x N identity matrix. The rank of this Hessian is N.

The Hessian of the network's function f(x) is H_f = sum_{j=1 to H} c_j * GeLU''(w_j^T * x + b_j) * (w_j * w_j^T).
Each term w_j * w_j^T is a rank-1 matrix. The sum of H such matrices can have a rank of at most H.
For the network to approximate g(x) with arbitrary precision, its Hessian H_f must be able to approximate H_g. This implies that the rank of H_f must be at least the rank of H_g.
Therefore, we must have H >= rank(H_g) = N. This gives us a lower bound.

Step 2: Analyzing the 1D case (N=1)

For N=1, we want to approximate g(x) = x^2. The lower bound from Step 1 is H >= 1.
Let's check if H=1 is sufficient. Can f(x) = c * GeLU(w*x + b) approximate x^2?
The function x^2 is symmetric around x=0. The GeLU function is not symmetric. A scaled and shifted GeLU function cannot be made symmetric to approximate x^2 over any significant domain.
Furthermore, the second derivative of g(x) is g''(x) = 2 (a constant), while the second derivative of f(x) is f''(x) = c*w^2*GeLU''(w*x+b), which is not constant.
So, H=1 is not sufficient.

To create a symmetric function, we can combine two GeLU neurons. Consider the function:
f(x) = c * (GeLU(a*x) + GeLU(-a*x))
This function is even (symmetric), since f(-x) = f(x). Using the Taylor expansion of GeLU around 0, we can show that for small x, this function is approximately proportional to x^2. This construction uses H=2 neurons.
So for N=1, the minimum width is 2.

Step 3: A constructive proof for N dimensions

We want to approximate g(x) = sum_{i=1 to N} x_i^2.
Based on the 1D case, we can approximate each x_i^2 term using two neurons.
For each dimension i in {1, ..., N}, we introduce two hidden neurons:
1. Neuron 1: Has weights w = a*e_i (where e_i is the standard basis vector for the i-th dimension and 'a' is a scaling parameter) and zero bias. It computes GeLU(a*x_i).
2. Neuron 2: Has weights w = -a*e_i and zero bias. It computes GeLU(-a*x_i).

The output layer can then compute a weighted sum of the activations from these neurons.
By choosing the output weights appropriately, we can form the sum:
f(x) = sum_{i=1 to N} c * (GeLU(a*x_i) + GeLU(-a*x_i))
This function approximates:
f(x) approx sum_{i=1 to N} C * x_i^2 = C * (sum_{i=1 to N} x_i^2)
where C is some constant. This provides a valid approximation of the squared norm.

This construction requires 2 neurons for each of the N dimensions. The weight vectors for these neurons (+-a*e_i) are all distinct. Thus, the total number of hidden neurons required is H = 2 * N.

Step 4: Conclusion

From Step 1, we have a lower bound H >= N.
From Step 2, we know that for N=1, H=2, which matches the formula 2*N.
From Step 3, we have a constructive proof that H = 2*N neurons are sufficient.
The lower bound H>=N can be tightened by observing that H=N is not sufficient, because the Hessian of the network cannot be made constant. The symmetry argument for requiring two neurons per squared term is the most direct path to the solution.

Therefore, the minimum hidden-layer width required is 2N.
"""
    print(explanation)
    # The prompt requires printing the equation "numbers"
    print("The final equation is: H = 2 * N")


if __name__ == '__main__':
    solve_neural_network_width()