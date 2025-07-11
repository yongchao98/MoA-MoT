import numpy as np

# This code block demonstrates the logic of the argument for the Heavy-ball method.
# It shows that if convergence to a point x_star is assumed, the limit point must be stationary.
# The existence of counter-examples relies on subtle behaviors not captured by this simple limit argument,
# but this illustrates the standard analysis.

# Let's define the terms in the rearranged heavy-ball equation:
# gamma * grad_f(x_k) = (x_k - x_{k+1}) + beta * (x_k - x_{k-1})

# Assume the sequence x_k converges to a point x_star.
# For demonstration, let's create a sequence that converges, e.g., x_k = 1 + (1/2)^k
# So, x_star = 1.
k = np.arange(10, 20) # Use large k to simulate k -> infinity
x_k = 1 + (0.5)**k
x_k_minus_1 = 1 + (0.5)**(k-1)
x_k_plus_1 = 1 + (0.5)**(k+1)

# Let's pick a value for beta
beta = 0.9
gamma = 0.1

# Calculate the two terms on the right-hand side (RHS) of the equation
term1 = x_k - x_k_plus_1
term2 = beta * (x_k - x_k_minus_1)

# Calculate the full RHS
RHS = term1 + term2

# According to the equation, this must be equal to the left-hand side (LHS)
# LHS = gamma * grad_f(x_k)
# So, grad_f(x_k) = RHS / gamma
grad_f_x_k = RHS / gamma

# Now, let's print the values as k gets large to see what they converge to
print("Demonstration of the limit argument for the Heavy-ball method:")
print("Let's assume x_k converges to x_star = 1, following the sequence x_k = 1 + 0.5^k")
print("-" * 20)
print(f"{'k':<5}{'x_k':<15}{'grad_f(x_k)':<15}")
for i in range(len(k)):
    print(f"{k[i]:<5}{x_k[i]:<15.10f}{grad_f_x_k[i]:<15.10f}")

print("-" * 20)
print("As k -> infinity, x_k -> 1.")
print("The term (x_k - x_{k+1}) -> 0.")
print("The term beta * (x_k - x_{k-1}) -> 0.")
print("Therefore, the right-hand side of the equation must converge to 0.")
print("This implies that gamma * grad_f(x_k) must also converge to 0.")
print("Since gamma > 0, grad_f(x_k) must converge to 0.")
print("If grad_f is continuous, grad_f(x_star) must be 0.")
print("Final equation if x_k converges to x_star:")
print(f"{gamma} * grad_f(x_star) = (x_star - x_star) + {beta} * (x_star - x_star)")
print(f"{gamma} * grad_f(x_star) = 0 + {beta} * 0")
print(f"{gamma} * grad_f(x_star) = 0")
print("grad_f(x_star) = 0")
print("\nConclusion: While this analysis shows convergence implies stationarity, the Heavy-ball method (3) is known to have pathological counterexamples, unlike (1) and (2).")