import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st

def gelu(x):
    """Gaussian Error Linear Unit."""
    return x * st.norm.cdf(x)

def gelu_prime(x):
    """First derivative of GeLU."""
    return st.norm.cdf(x) + x * st.norm.pdf(x)

def gelu_double_prime(x):
    """Second derivative of GeLU."""
    # d/dx(Φ(x) + xφ(x)) = φ(x) + (φ(x) + x(-xφ(x))) = 2φ(x) - x^2φ(x)
    # The above is incorrect, d/dx(xφ(x)) = φ(x) - x^2 φ(x)
    # The correct one is d/dx(x*phi(x)) = phi(x) + x*(-x*phi(x)) = phi(x) - x^2*phi(x)
    # So d/dx(phi(x) + x*phi(x)) = phi(x) + phi(x) - x^2*phi(x) = 2phi(x) - x^2phi(x). WRONG.
    # d/dx[x*st.norm.pdf(x)] = st.norm.pdf(x) + x * (-x*st.norm.pdf(x))
    return st.norm.pdf(x) + st.norm.pdf(x) + x * (-x * st.norm.pdf(x)) # This is also wrong
    
    # Let's use the verified formula: GeLU''(x) = φ(x) + (φ(x) - x^2φ(x)) is not right.
    # GeLU''(x) = φ(x) + (d/dx (xφ(x))) = φ(x) + (φ(x) + x*φ'(x)) = 2φ(x) -x^2 φ(x) also wrong
    # Let's re-re-derive:
    # d/dx [Φ(x) + xφ(x)] = φ(x) + [1*φ(x) + x*(-x*φ(x))] = 2φ(x) - x^2φ(x) -- No, this is what I had before. Wait.
    # My calculus is rusty.
    # d/dx (x * φ(x)) = 1*φ(x) + x*φ'(x) = φ(x) + x*(-x*φ(x)) = φ(x)*(1-x^2)
    # So, GeLU''(x) = d/dx[Φ(x)] + d/dx[xφ(x)] = φ(x) + φ(x)*(1-x^2) = φ(x)*(2-x^2).
    # Okay, this must be it.
    return st.norm.pdf(x) * (2 - x**2)


def approx_square(z, b_offset, c1, c2):
    """Approximates z^2 using 2 GeLU neurons."""
    # The approximation is c1 * [GeLU(z+b) + GeLU(-z+b)] + c2
    return c1 * (gelu(z + b_offset) + gelu(-z + b_offset)) + c2

# --- Parameters for approximation ---
# Based on Taylor expansion: z^2 ≈ (1/GeLU''(b)) * [h(z) - 2GeLU(b)]
# So c1 = 1/GeLU''(b) and c2 = -2*GeLU(b)*c1
b = 1.5  # A chosen offset
c1 = 1.0 / gelu_double_prime(b)
c2 = -2 * gelu(b) * c1

# Generate input data
z = np.linspace(-2, 2, 400)
# Calculate the true squared values
true_square = z**2
# Calculate the approximated squared values
approximated_square = approx_square(z, b, c1, c2)

# --- Output ---
# The core result is that the minimum hidden-layer width is 2 * N.
# We demonstrate how 2 neurons can approximate z^2, which is the building block.
print("To compute the squared norm of an N-dimensional vector, a shallow neural network with GeLU activations requires a hidden layer of width 2N.")
print("This is because the function ||x||^2 = x_1^2 + ... + x_N^2 can be decomposed into N separate quadratic problems.")
print("Each quadratic problem, z^2, can be approximated to arbitrary precision using 2 GeLU neurons.")
print("The construction for z^2 is based on the approximation: z^2 ≈ c1 * (GeLU(z+b) + GeLU(-z+b)) + c2")
print("\nDemonstration for N=1 (approximating z^2):")
print(f"Using parameters: b={b:.2f}, c1={c1:.2f}, c2={c2:.2f}")
print("Approximation formula: y ≈ {:.2f} * (GeLU(z + {:.2f}) + GeLU(-z + {:.2f})) + ({:.2f})".format(c1, b, b, c2))
print("Therefore, combining the approximations for each dimension x_i requires N * 2 = 2N neurons.")

# Plotting the result for visualization
plt.figure(figsize=(8, 6))
plt.plot(z, true_square, label='True function: $z^2$', linewidth=2)
plt.plot(z, approximated_square, label=f'Approximation with 2 GeLU neurons', linestyle='--', linewidth=2)
plt.title('Approximating $z^2$ with 2 GeLU Neurons')
plt.xlabel('z')
plt.ylabel('Value')
plt.legend()
plt.grid(True)
plt.show()
