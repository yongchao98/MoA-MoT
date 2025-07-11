import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.special import erfc

def plot_smooth_V_shape():
    """
    This function demonstrates that statement B is true by constructing and plotting
    a smooth curve gamma(t) whose image is the set L = {(x, y) | y = |x|}.
    """

    # We construct a smooth bijection x: R -> R with a zero of infinite order at t=0.
    # We use the scaled and shifted complementary error function, which is smooth and
    # provides a nice transition from 0 to 1.
    # Let phi(t) = 0.5 * erfc(-t). This is a smooth function from R to (0, 1).
    # To get a map from R to R, we can use the tan function on a scaled version.
    # Let h(t) = tan(pi * (phi(t) - 0.5)). This maps R to R smoothly.
    # x(t) = h(t) will be a smooth bijection R -> R.
    # Let's check the derivatives at t=0. The derivatives of erfc(-t) at t=0 are non-zero.
    # This construction is a bit tricky.
    
    # A simpler construction for visualization uses a function known to have a zero of infinite order:
    # f(t) = exp(-1/t^2) for t!=0, and f(0)=0.
    # To get a bijection R->R, we integrate this function and then use tan.
    
    # Define the integrand for our smooth bijection
    def integrand(s):
        if s == 0:
            return 0
        return np.exp(-1 / s**2)

    # Calculate the definite integral from -inf to inf, to use as a normalization constant
    c, _ = quad(integrand, -np.inf, np.inf)

    # Define g(t) = (pi/c) * integral from 0 to t of the integrand.
    # This gives a smooth map from R to (-pi/2, pi/2).
    def g(t):
        val, _ = quad(integrand, 0, t)
        return (np.pi / c) * val

    # Vectorize g for plotting
    g_vec = np.vectorize(g)

    # Define x(t) = tan(g(t)). This is a smooth bijection from R to R.
    # All derivatives of x(t) at t=0 are zero.
    def x_func(t):
        return np.tan(g_vec(t))

    # The curve is gamma(t) = (x(t), |x(t)|)
    t = np.linspace(-4, 4, 500)
    # The t=0 point can cause numerical instability, let's handle it
    t = t[t != 0]

    x_vals = x_func(t)
    y_vals = np.abs(x_vals)

    # Create the plot
    plt.figure(figsize=(8, 8))
    plt.plot(x_vals, y_vals, label=r'Smooth curve $\gamma(t)$')
    # For comparison, plot the actual y=|x|
    x_comp = np.linspace(min(x_vals), max(x_vals), 100)
    y_comp = np.abs(x_comp)
    plt.plot(x_comp, y_comp, 'r--', label=r'$y=|x|$ (for reference)')
    
    plt.title(r'Visualization of a smooth curve $\gamma(\mathbb{R})$ tracing $L = \{(x,y) : y=|x|\}$')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.axhline(0, color='black',linewidth=0.5)
    plt.axvline(0, color='black',linewidth=0.5)
    plt.grid(True, linestyle=':')
    plt.axis('equal')
    plt.legend()
    
    print("This plot shows a smooth curve whose image is the set y=|x|, demonstrating that statement B is true.")
    print("The false statement is C.")
    plt.show()

# Run the function to generate the plot
plot_smooth_V_shape()
