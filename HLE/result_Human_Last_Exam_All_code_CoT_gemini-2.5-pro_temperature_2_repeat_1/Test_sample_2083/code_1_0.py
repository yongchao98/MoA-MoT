import math

def calculate_minimum_width(N):
  """
  Calculates the minimum hidden-layer width for a shallow neural network
  to compute the squared norm of an N-dimensional vector.

  The target function is f(x) = ||x||^2 = x_1^2 + x_2^2 + ... + x_N^2.

  The reasoning is as follows:
  1. The function f(x) is a sum of N independent quadratic functions.
  2. Each 1D quadratic function, g(z) = z^2, is not monotonic. A single GELU
     neuron produces a monotonic output, so at least 2 neurons are needed
     to approximate z^2.
  3. A combination of two neurons, GELU(w*z) + GELU(-w*z), creates an even
     function that can approximate z^2 on a compact domain.
  4. To compute the full function sum(x_i^2), we need to approximate each x_i^2
     term independently. Using shared neurons would introduce unwanted cross-terms
     (e.g., x_i * x_j).
  5. Therefore, we need 2 neurons for each of the N dimensions.
  6. The minimum required width is 2 * N.
  """

  # The minimum hidden-layer width is 2 * N.
  width = 2 * N
  
  print(f"For an input vector of dimension N = {N}, the minimum hidden-layer width is 2 * N.")
  print(f"Final Equation: 2 * {N} = {width}")

# Example with N = 100
n_dimensions = 100
calculate_minimum_width(n_dimensions)
