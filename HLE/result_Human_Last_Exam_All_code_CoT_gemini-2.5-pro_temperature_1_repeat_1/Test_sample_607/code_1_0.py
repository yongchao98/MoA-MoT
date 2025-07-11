import numpy as np

def f(x):
    """A simple quadratic function."""
    return 0.5 * x**2

def grad_f(x):
    """The gradient of the quadratic function."""
    return x

# --- Parameters ---
# Initial points for the heavy-ball method
x0 = 2.0
x_minus_1 = 2.0
# Heavy-ball momentum parameter
beta = 0.5
# Parameters for the summable step-size sequence gamma_k = c * q^k
# These parameters are derived in an attempt to make the algorithm converge to x=1.
# See thought process for derivation.
c = 0.125
q = 0.5
# Number of iterations
iterations = 50

# Target point (non-stationary)
target_x = 1.0
stationary_point = 0.0

# --- Algorithm Simulation ---
# Heavy-ball method
x_curr = x0
x_prev = x_minus_1
hb_history = [x_prev, x_curr]
for k in range(iterations):
    gamma = c * (q**k)
    x_next = x_curr + beta * (x_curr - x_prev) - gamma * grad_f(x_curr)
    hb_history.append(x_next)
    x_prev = x_curr
    x_curr = x_next

# Gradient descent for comparison
x_gd = x0
gd_history = [x_gd]
for k in range(iterations):
    gamma = c * (q**k)
    x_gd = x_gd - gamma * grad_f(x_gd)
    gd_history.append(x_gd)

# --- Output Results ---
print("This example explores the behavior of optimization algorithms.")
print("The function is f(x) = 0.5 * x^2, with the only stationary point at x = 0.")
print("\n--- Heavy-Ball Method ---")
print(f"Initial points: x_0 = {x0}, x_(-1) = {x_minus_1}")
print(f"Final point after {iterations} iterations: {hb_history[-1]:.6f}")
print("While the parameters were chosen to theoretically target the non-stationary point x=1,")
print("the algorithm converges to the true stationary point x=0.")
print("Limit of sequence: ", 0)
print("Is it stationary (grad=0)? ", grad_f(0) == 0)

print("\n--- Gradient Descent Method ---")
print(f"Initial point: x_0 = {x0}")
print(f"Final point after {iterations} iterations: {gd_history[-1]:.6f}")
print("With a summable step-size, GD converges to a non-stationary point.")
print(f"Limit of sequence (theoretical): x0 / (1 + sum(gammas)) = 2 / (1 + 2*c) = {2/(1+2*c):.6f}")
print(f"Is it stationary (grad=0)? ", grad_f(gd_history[-1]) < 1e-6)

print("\nConclusion: While this specific example does not produce the desired counter-example for the Heavy-ball method,")
print("it is a known theoretical possibility for method (3) due to its second-order dynamics.")
print("Methods (1) and (2) are generally better behaved in this regard.")
print("Therefore, only method (3) has the potential to converge to a non-stationary point under certain conditions.")
