import sympy

def demonstrate_moments():
    """
    This function demonstrates that for a non-zero Schwartz function,
    like the Gaussian f(x) = exp(-x**2), the moments are not all zero.
    This supports the theorem that if all moments are zero, the function
    must be the zero function.
    """
    # Define the symbolic variable and the function
    x = sympy.Symbol('x')
    f = sympy.exp(-x**2)

    print(f"Demonstration using the non-zero Schwartz function: f(x) = {f}")
    print("The k-th moment is M_k = integral(x**k * f(x) dx) from -oo to +oo.")
    print("-" * 50)

    # Calculate and print moments for k = 0 to 7
    for k in range(8):
        # The term inside the integral
        integrand = x**k * f

        # Perform the definite symbolic integration from -oo to +oo
        moment_value = sympy.integrate(integrand, (x, -sympy.oo, sympy.oo))

        # Output the equation for the calculated moment
        # We output each number involved in the calculation as requested
        print(f"Equation for k={k}: M_{k} = Integral(x**{k} * exp(-x**2), x=-oo..oo) = {moment_value}")

    print("-" * 50)
    print("As shown, the moments for even k are non-zero.")
    print("This is consistent with the theorem: a non-zero function must have non-zero moments.")

if __name__ == '__main__':
    demonstrate_moments()