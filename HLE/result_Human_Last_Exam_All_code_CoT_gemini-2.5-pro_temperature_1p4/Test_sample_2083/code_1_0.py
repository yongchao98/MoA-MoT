import sympy

# Define N as a symbolic variable
N = sympy.Symbol('N')

# According to the step-by-step derivation:
# 1. The target function is the sum of N independent quadratic terms (x_j^2).
# 2. The GeLU activation function is non-monotonic and U-shaped, just like x^2.
# 3. Therefore, one GeLU neuron is sufficient to approximate each x_j^2 term.
# 4. By dedicating one neuron to each of the N dimensions, we can construct the sum.
# 5. The minimum required hidden-layer width is N.

min_hidden_width = N

# The final equation is simply H = N.
# We will print the variables in the equation as requested.
print(f"The minimum hidden-layer width (H) is a function of the input dimension (N).")
print(f"The required width H is equal to N.")
print(f"Final Equation: H = {min_hidden_width}")
