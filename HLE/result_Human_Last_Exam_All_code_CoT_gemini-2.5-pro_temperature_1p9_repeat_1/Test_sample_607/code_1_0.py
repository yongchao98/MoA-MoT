import numpy as np

# We choose the Heavy-ball method (3) to demonstrate the possibility of
# converging to a non-stationary point.
# The strategy is to construct a sequence of points {x_k} that converges to a
# chosen non-stationary point, and then from the algorithm's update rule,
# we solve for the step sizes {gamma_k} that must have been used. If these
# step sizes are valid (i.e., positive), then we have shown that it's possible.

# 1. Define a smooth function and a target non-stationary point.
def f(x):
    """A simple quadratic function."""
    return 0.5 * x**2

def nabla_f(x):
    """The gradient of f(x)."""
    return x

# The unique stationary point is x=0, where nabla_f(x) = 0.
# We will target convergence to a non-stationary point, x_star = 1.
x_star = 1.0

# 2. Define algorithm parameters.
beta = 0.5  # A standard choice for the momentum parameter.

# 3. Define the desired sequence x_k that converges to x_star.
# We can use a geometric progression for the error term.
# x_k = x_star + c * rho^k, with |rho| < 1 for convergence.
c = 1.0
rho = 0.5
num_steps = 30

# Generate the ideal sequence from k=-1 to k=num_steps. We need x_{k-1} to start.
# This sequence converges to x_star = 1.
x_sequence_ideal = np.array([x_star + c * rho**k for k in range(-1, num_steps + 1)])


# 4. Reverse-engineer the step sizes gamma_k.
# From x_{k+1} = x_k + beta*(x_k - x_{k-1}) - gamma_k * nabla_f(x_k), we get:
# gamma_k = (x_k + beta*(x_k - x_{k-1}) - x_{k+1}) / nabla_f(x_k)
gamma_sequence = []
for k in range(num_steps):
    # For iteration k, we use x_k, x_{k-1} from our ideal sequence
    # to compute x_{k+1}, also from the ideal sequence.
    # The array indices are shifted: x_k is at index k+1, x_{k-1} at k, x_{k+1} at k+2.
    x_k = x_sequence_ideal[k+1]
    x_k_minus_1 = x_sequence_ideal[k]
    x_k_plus_1 = x_sequence_ideal[k+2]

    numerator = x_k + beta * (x_k - x_k_minus_1) - x_k_plus_1
    denominator = nabla_f(x_k)

    # For x_star=1, nabla_f(x_k)=x_k is always positive and far from 0.
    gamma_k = numerator / denominator
    gamma_sequence.append(gamma_k)

# 5. Run a simulation with these step sizes to prove the concept.
print("--- Heavy-ball method converging to a non-stationary point ---")
print(f"Function f(x) = 0.5 * x^2, Stationary point at x=0.")
print(f"Targeting convergence to non-stationary point x* = {x_star}")
print(f"Gradient at target: nabla_f({x_star}) = {nabla_f(x_star)}\n")

# Initialize the algorithm with the first two points of our ideal sequence.
x_previous = x_sequence_ideal[0]
x_current = x_sequence_ideal[1]

print(f"k=0,  x_k = {x_current:.6f}")

for k in range(num_steps - 1):
    # Use the pre-calculated step size for this iteration
    gamma_k = gamma_sequence[k]
    
    # Apply the Heavy-ball update
    x_next = x_current + beta * (x_current - x_previous) - gamma_k * nabla_f(x_current)
    
    x_previous = x_current
    x_current = x_next
    
    if (k + 1) % 5 == 0:
        print(f"k={k+1:<2}, x_k = {x_current:.6f}, (gamma_k = {gamma_k:.6f})")

print(f"\nAfter {num_steps-1} iterations, the algorithm is at x = {x_current:.6f}.")
print(f"This demonstrates convergence towards the non-stationary point {x_star}.")
