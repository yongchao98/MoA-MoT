import numpy as np

# This script demonstrates that Gradient Descent (1) can converge to a 
# non-stationary point if a summable step size sequence is used.
# The same principle applies to algorithms (2) and (3).

# Define a simple smooth function and its gradient.
# Let f(x) = 0.5 * x^2.
# The unique stationary point is x = 0, where the gradient is zero.
def f(x):
    """The function to minimize."""
    return 0.5 * x**2

def grad_f(x):
    """The gradient of the function."""
    return x

# --- Algorithm Parameters ---
# Initial point
x0 = 10.0
# Number of iterations
num_iterations = 100
# Parameters for the summable step size sequence: gamma_k = a * q^k
# The sum of this geometric series is finite, so the steps die out quickly.
a = 0.5
q = 0.5

# --- Main Logic ---
print("--- Demonstration of Convergence to a Non-Stationary Point ---")
print(f"Function f(x) = 0.5 * x^2")
print(f"Gradient grad(f)(x) = x")
print(f"The only stationary point for this function is x = 0.")
print(f"\nStarting at initial point x_0 = {x0}")
print(f"Using a summable step size sequence gamma_k = {a} * ({q})^k.")

# Initialize the iterate
x = x0

# Run the Gradient Descent algorithm
for k in range(num_iterations):
    # Calculate the current step size
    gamma_k = a * (q**k)
    # Calculate the gradient at the current point
    gradient = grad_f(x)
    
    # Apply the update rule: x_{k+1} = x_k - gamma_k * grad(x_k)
    x_new = x - gamma_k * gradient
    
    # As requested, output the numbers for the first step's equation
    if k == 0:
        print("\n--- Example: First Iteration (k=0) ---")
        print(f"Equation: x_1 = x_0 - gamma_0 * grad(x_0)")
        print(f"Numbers:  x_1 = {x:.4f} - {gamma_k:.4f} * {gradient:.4f}")
        print(f"Result:   x_1 = {x_new:.4f}\n")
        
    # Update the iterate
    x = x_new

# The algorithm converges after many iterations
final_x = x
final_grad = grad_f(final_x)

# --- Final Result and Conclusion ---
print(f"--- Final Result after {num_iterations} iterations ---")
# The theoretical limit is non-zero because x_k = x_0 * Product(1-gamma_i)
# and since Sum(gamma_i) is finite, the product converges to a non-zero number.
print(f"The sequence converges to the point x* = {final_x:.8f}")

# Check for stationarity at the limit point
print("\nChecking if the limit point x* is stationary:")
print(f"The gradient at x* is grad_f(x*) = {final_grad:.8f}")

# A point is stationary if its gradient is zero.
# We use a small tolerance for floating-point comparison.
if np.isclose(final_grad, 0.0):
    print("\nConclusion: The limit point IS stationary.")
else:
    print("\nConclusion: The limit point IS NOT stationary.")
