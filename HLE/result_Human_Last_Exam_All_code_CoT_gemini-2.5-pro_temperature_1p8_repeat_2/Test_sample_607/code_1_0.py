import numpy as np

def rosenbrock(x):
    """The Rosenbrock function."""
    return (1 - x[0])**2 + 100 * (x[1] - x[0]**2)**2

def rosenbrock_grad(x):
    """The gradient of the Rosenbrock function."""
    grad = np.zeros_like(x)
    grad[0] = -2 * (1 - x[0]) - 400 * x[0] * (x[1] - x[0]**2)
    grad[1] = 200 * (x[1] - x[0]**2)
    return grad

def heavy_ball_method(grad_f, x0, x_prev, gamma, beta, max_iter=1000, tol=1e-6):
    """
    Heavy-ball method for optimization.
    
    x_k+1 = x_k + beta * (x_k - x_k-1) - gamma * grad_f(x_k)
    """
    x_curr = x0
    path = [x_curr]
    
    print("Running Heavy-ball method...")
    for k in range(max_iter):
        grad = grad_f(x_curr)
        x_next = x_curr + beta * (x_curr - x_prev) - gamma * grad
        
        path.append(x_next)
        
        if np.linalg.norm(x_next - x_curr) < tol:
            print(f"Converged after {k+1} iterations.")
            break
            
        x_prev = x_curr
        x_curr = x_next
        
    return x_curr, np.array(path)

# Parameters for the Heavy-ball method
x_start = np.array([-1.5, -1.0])
# We need two initial points for the heavy-ball method
x0 = x_start
x_prev = x_start 

# Learning rate and momentum parameter
gamma = 0.001 
beta = 0.9    

# Run the optimization
final_x, path = heavy_ball_method(rosenbrock_grad, x0, x_prev, gamma, beta)

# The stationary point for Rosenbrock function is at (1, 1)
# The gradient at this point is [0, 0]
stationary_point = np.array([1.0, 1.0])

print(f"\nAlgorithm (3) - Heavy-ball Method:")
print(f"Starting point: {x_start}")
print(f"Final point reached: [{final_x[0]:.6f}, {final_x[1]:.6f}]")
grad_at_final = rosenbrock_grad(final_x)
print(f"Gradient at final point: [{grad_at_final[0]:.6f}, {grad_at_final[1]:.6f}]")
print(f"Norm of gradient at final point: {np.linalg.norm(grad_at_final):.6f}")

print("\nConclusion:")
print("Although a straightforward analysis suggests all three methods converge to stationary points,")
print("the Heavy-ball method (3) is known from optimization literature to have counter-examples")
print("where it can converge to a non-stationary point for specific non-convex functions.")
print("The other two methods are generally considered safe in this regard under standard assumptions.")
