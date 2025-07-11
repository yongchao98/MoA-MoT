import sympy as sp

def demonstrate_moment_fourier_relation():
    """
    This function demonstrates the mathematical relationship at the heart of the proof.
    It shows that the k-th derivative of the Fourier Transform of a function f(x)
    at xi=0 is proportional to the k-th moment of f(x). If all moments are zero,
    then all derivatives of the Fourier transform at the origin are zero.
    """
    # Define symbolic variables
    x, xi = sp.symbols('x xi', real=True)
    f_gaussian = sp.exp(-sp.pi * x**2)

    # In Sympy, the Fourier transform of exp(-pi*x**2) is exp(-pi*xi**2)
    F_gaussian = sp.fourier_transform(f_gaussian, x, xi, noconds=True)

    print("We consider the Schwartz function f(x) = exp(-pi*x**2).")
    print(f"Its Fourier transform is F(xi) = {F_gaussian}\n")

    print("The proof relies on the identity:")
    print("d^k F / d(xi)^k |_(xi=0) = (-2*pi*I)**k * Integral(x**k * f(x) dx)")
    print("Let's verify this for k = 0, 1, 2, 3, 4.\n")

    # The condition is that all moments are zero. Here we calculate them for our example function.
    # If we were to assume they are zero, the RHS would be zero, forcing the LHS to be zero.
    for k_val in range(5):
        # Calculate the k-th moment M_k of f(x)
        moment_k = sp.integrate(x**k_val * f_gaussian, (x, -sp.oo, sp.oo))

        # Calculate the k-th derivative of the Fourier transform F(xi)
        F_deriv_k = sp.diff(F_gaussian, xi, k_val)

        # Evaluate the derivative at xi = 0
        F_deriv_k_at_0 = F_deriv_k.subs(xi, 0)

        # The proportionality constant from the formula
        C_k = (-2 * sp.pi * sp.I)**k_val
        
        # Calculate the right-hand side of the identity
        rhs = sp.simplify(C_k * moment_k)

        print(f"--- For k = {k_val} ---")
        print(f"Moment M_{k_val} = {moment_k}")
        print(f"LHS: d^{k_val}F/d(xi)^{k_val} at xi=0 = {F_deriv_k_at_0}")
        print(f"RHS: (-2*pi*I)**{k_val} * M_{k_val} = {rhs}")
        
        # Verify the equality
        is_equal = sp.simplify(F_deriv_k_at_0 - rhs) == 0
        print(f"Are LHS and RHS equal?  {is_equal}")
        print("-" * 30)

if __name__ == '__main__':
    demonstrate_moment_fourier_relation()