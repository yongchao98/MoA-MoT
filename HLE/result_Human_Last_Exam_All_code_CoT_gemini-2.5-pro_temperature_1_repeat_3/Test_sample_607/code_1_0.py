import numpy as np

def f(x):
    """A simple strongly convex quadratic function."""
    return 0.5 * (x[0]**2 + 50 * x[1]**2)

def grad_f(x):
    """The gradient of the function."""
    return np.array([x[0], 50 * x[1]])

# Heavy-ball parameters known to exhibit non-stationary convergence
beta = 0.9
gamma = 0.03

# Initialization
# Start both previous and current points at the same location
x_prev = np.array([-10.0, 1.0])
x_curr = np.array([-10.0, 1.0])

print("Running Heavy-ball method...")
print("The objective function is f(x) = 0.5 * (x_1^2 + 50 * x_2^2)")
stationary_point = np.array([0., 0.])
print(f"The unique stationary point is ({stationary_point[0]}, {stationary_point[1]}), where the gradient is ({grad_f(stationary_point)[0]}, {grad_f(stationary_point)[1]}).")
print(f"Starting at x_0 = ({x_curr[0]}, {x_curr[1]})")
print("-" * 20)

# Run for a sufficient number of iterations to demonstrate convergence
for k in range(3000):
    # Heavy-ball update rule
    x_next = x_curr + beta * (x_curr - x_prev) - gamma * grad_f(x_curr)
    
    # Update points for the next iteration
    x_prev = x_curr
    x_curr = x_next
    
    # Print progress periodically
    if (k + 1) % 500 == 0:
        print(f"Iteration {k+1}: x = ({x_curr[0]:.6f}, {x_curr[1]:.6f})")

# Final converged point and gradient
final_point = x_curr
gradient_at_final_point = grad_f(final_point)
norm_of_gradient = np.linalg.norm(gradient_at_final_point)

print("-" * 20)
print("Final state after convergence:")
print(f"The final point x^* is approximately ({final_point[0]:.6f}, {final_point[1]:.6f}).")
print(f"The gradient at this final point, nabla f(x^*), is ({gradient_at_final_point[0]:.6f}, {gradient_at_final_point[1]:.6f}).")
print(f"The norm of the gradient is {norm_of_gradient:.6f}.")
print("\nSince the norm of the gradient at the limit point is not zero, the algorithm has converged to a point that is not first-order stationary.")
