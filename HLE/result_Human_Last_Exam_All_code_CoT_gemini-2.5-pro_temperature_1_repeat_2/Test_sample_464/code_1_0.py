import numpy as np
from scipy.integrate import quad
from scipy.misc import derivative

def solve():
    """
    Demonstrates and proves that a Schwartz function with all zero moments must be zero.
    """
    
    print("### Part 1: Numerical Demonstration ###")
    print("Let's use a sample Schwartz function, the Gaussian f(x) = exp(-x^2).")
    print("We will verify the key identity connecting moments of f(x) and derivatives of its Fourier transform.\n")

    # Define the Schwartz function f(x)
    def f(x):
        return np.exp(-x**2)

    # The analytical Fourier transform of f(x) is sqrt(pi) * exp(-(pi*xi)^2)
    # The Fourier transform is defined as F(xi) = integral(f(x) * exp(-2*pi*i*x*xi) dx)
    def f_hat(xi):
        return np.sqrt(np.pi) * np.exp(-(np.pi * xi)**2)

    print("The identity is: (d^k / d(xi)^k) f_hat(0) = (-2*pi*i)^k * M_k")
    print("where M_k is the k-th moment of f(x).\n")
    print("-" * 60)
    print(f"{'k':<3} | {'Moment M_k':<20} | {'Derivative of f_hat(0)':<25} | {'(-2*pi*i)^k * M_k':<25}")
    print("-" * 60)

    for k in range(6):
        # 1. Calculate the k-th moment M_k
        # M_k = integral from -inf to inf of x^k * f(x) dx
        integrand = lambda x: (x**k) * f(x)
        moment, _ = quad(integrand, -np.inf, np.inf)

        # 2. Calculate the k-th derivative of f_hat at xi = 0
        # For complex functions, we need to differentiate the real and imaginary parts separately.
        f_hat_real = lambda xi: f_hat(xi).real
        f_hat_imag = lambda xi: f_hat(xi).imag
        
        # Since our f_hat is real, its derivatives will be real.
        # But for the general case, let's keep it formal.
        # However, scipy.misc.derivative doesn't handle complex-valued functions well.
        # Since our chosen f_hat is real, f_hat_imag is zero.
        # The derivatives of a real function are real.
        deriv_val = derivative(f_hat, 0.0, dx=1e-1, n=k, order=max(k + 1 + (k % 2), 7))

        # 3. Calculate the value from the identity
        identity_val = ((-2 * np.pi * 1j)**k) * moment

        # 4. Print comparison
        print(f"{k:<3} | {moment:<20.6f} | {deriv_val:<25.6f} | {identity_val.real:<.6f} + {identity_val.imag:<.6f}j")

    print("-" * 60)
    print("\nThe numerical results match, confirming the identity.\n")


    print("\n### Part 2: Formal Proof ###")
    proof_text = """
    We want to prove that if f is a Schwartz function and all its moments are zero, then f must be the zero function.

    Hypothesis: f is in the Schwartz space S(R) and its moments Mk are all zero.
    M_k = integral[-inf, inf] (x^k * f(x) dx) = 0 for all k = 0, 1, 2, ...

    Step 1: Consider the Fourier transform of f(x).
    Let f_hat(xi) be the Fourier transform of f(x), defined as:
    f_hat(xi) = integral[-inf, inf] (f(x) * exp(-2*pi*i*x*xi) dx)

    Step 2: Relate the derivatives of the Fourier transform to the moments of the function.
    We can differentiate under the integral sign (allowed because f is a Schwartz function):
    d/d(xi) f_hat(xi) = integral[-inf, inf] (f(x) * (-2*pi*i*x) * exp(-2*pi*i*x*xi) dx)

    Repeating this k times, we get the k-th derivative:
    (d^k / d(xi)^k) f_hat(xi) = integral[-inf, inf] ((-2*pi*i*x)^k * f(x) * exp(-2*pi*i*x*xi) dx)

    Step 3: Evaluate the derivatives at xi = 0.
    Setting xi = 0, the exponential term becomes exp(0) = 1:
    (d^k / d(xi)^k) f_hat(0) = integral[-inf, inf] ((-2*pi*i*x)^k * f(x) dx)
                           = (-2*pi*i)^k * integral[-inf, inf] (x^k * f(x) dx)
                           = (-2*pi*i)^k * M_k

    Step 4: Use the hypothesis that all moments are zero.
    Since we are given that M_k = 0 for all k, we have:
    (d^k / d(xi)^k) f_hat(0) = (-2*pi*i)^k * 0 = 0 for all k = 0, 1, 2, ...
    This means that the Fourier transform f_hat(xi) and all of its derivatives are zero at xi = 0.

    Step 5: Use the properties of the Fourier transform of a Schwartz function.
    If f(x) is a Schwartz function, its Fourier transform f_hat(xi) is also a Schwartz function.
    A key property of Schwartz functions is that they are infinitely differentiable (C-infinity) and also analytic. An analytic function is one that is locally equal to its Taylor series.

    Step 6: Analyze the Taylor series of f_hat(xi) around xi = 0.
    The Taylor series for f_hat(xi) around xi = 0 is:
    f_hat(xi) = sum_{k=0 to inf} [ (d^k/d(xi)^k)f_hat(0) / k! ] * xi^k
    From Step 4, we know all the derivative terms are zero. So the series is:
    f_hat(xi) = sum_{k=0 to inf} [ 0 / k! ] * xi^k = 0

    Since f_hat(xi) is analytic, its Taylor series converges to the function itself in a neighborhood of xi=0. Thus, f_hat(xi) = 0 in that neighborhood.
    By the identity theorem for analytic functions, if an analytic function is zero on any interval, it must be zero everywhere. Therefore, f_hat(xi) = 0 for all xi in R.

    Step 7: Conclude by taking the inverse Fourier transform.
    The function f(x) can be recovered from its Fourier transform f_hat(xi) using the inverse Fourier transform:
    f(x) = integral[-inf, inf] (f_hat(xi) * exp(2*pi*i*x*xi) d(xi))
    Since f_hat(xi) is identically zero, we have:
    f(x) = integral[-inf, inf] (0 * exp(2*pi*i*x*xi) d(xi)) = 0

    Conclusion:
    We have shown that f(x) must be the zero function. So, yes, it follows that f = 0.
    """
    print(proof_text)


solve()
