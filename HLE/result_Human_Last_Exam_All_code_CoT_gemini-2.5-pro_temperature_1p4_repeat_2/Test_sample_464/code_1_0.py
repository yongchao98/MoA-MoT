import sympy
from sympy import Function, Symbol, Integral, I, pi, oo, Derivative, Eq

def prove_moment_problem_solution():
    """
    This function provides a step-by-step proof for the problem:
    If a Schwartz class function f(x) has all its moments equal to zero,
    does it follow that f = 0?
    
    The proof uses the properties of the Fourier transform.
    """
    # Introduction to the problem and the method.
    print("The problem asks whether a Schwartz function f(x) must be the zero function")
    print("if all its moments, Integral(x**k * f(x) dx), are zero for all non-negative integers k.")
    print("The answer is yes. We can prove this using the Fourier transform.\n")

    # Define symbolic variables for the demonstration
    x = Symbol('x', real=True)
    xi = Symbol('xi', real=True)
    k = Symbol('k', integer=True, nonnegative=True)
    f = Function('f')(x)
    f_hat_func = Function('\u00e2\u0178')(xi) # Unicode for f with a hat

    # --- Step 1: The Fourier Transform and its Derivatives ---
    print("--- Step 1: Relate the Fourier Transform to the Moments ---")
    print("Let \u00e2\u0178(xi) denote the Fourier transform of f(x).")
    print(f"By definition: \u00e2\u0178(xi) = Integral(f(x) * exp(-2*\u03c0*I*x*xi), (x, -oo, oo))\n")

    print("Since f(x) is a Schwartz function, we can differentiate under the integral sign.")
    print("The k-th derivative of \u00e2\u0178(xi) with respect to xi is:")
    kth_deriv_expr = Derivative(f_hat_func, (xi, k))
    kth_deriv_integral = Integral(f * (-2*pi*I*x)**k * sympy.exp(-2*pi*I*x*xi), (x, -oo, oo))
    print(Eq(kth_deriv_expr, kth_deriv_integral, evaluate=False), "\n")
    
    print("Evaluating the k-th derivative at xi = 0 (since exp(0) = 1):")
    moment_k = Integral(x**k * f, (x, -oo, oo))
    kth_deriv_at_0_final = ((-2*pi*I)**k) * moment_k
    print(Eq(kth_deriv_expr.subs(xi, 0), kth_deriv_at_0_final, evaluate=False), "\n")

    # --- Step 2: Apply the Given Condition ---
    print("--- Step 2: Apply the Condition that All Moments are Zero ---")
    print("The problem states that all moments of f(x) are zero. The k-th moment is:")
    print(f"  {moment_k} = 0  (for all k = 0, 1, 2, ...)\n")
    
    print("Substituting this into our equation for the derivatives:")
    print(f"  {kth_deriv_expr.subs(xi, 0)} = ((-2*\u03c0*I)**k) * 0 = 0\n")
    
    print("This implies that \u00e2\u0178(xi) and all its derivatives are zero at xi = 0.\n")
    
    # --- Step 3: Use Analyticity of the Fourier Transform ---
    print("--- Step 3: The Taylor Series and Conclusion for the Fourier Transform ---")
    print("A key property is that the Fourier transform of a Schwartz function is also")
    print("a Schwartz function, which means it is analytic (infinitely differentiable and")
    print("equal to its Taylor series in a neighborhood of any point).\n")
    
    print("The Taylor series expansion of \u00e2\u0178(xi) around xi = 0 is:")
    print("  \u00e2\u0178(xi) = \u03a3_{k=0 to \u221e} [\u00e2\u0178^(k)(0) / k!] * xi^k")
    print("Since we found that \u00e2\u0178^(k)(0) = 0 for all k, the Taylor series is identically zero:")
    print("  \u00e2\u0178(xi) = 0 + 0\u00b7xi + 0\u00b7xi\u00b2/2! + ... = 0\n")
    
    print("Since \u00e2\u0178(xi) is analytic, it must be zero everywhere. So, \u00e2\u0178(xi) = 0 for all xi.\n")

    # --- Step 4: Conclusion for f(x) ---
    print("--- Step 4: Final Conclusion for f(x) ---")
    print("The Fourier transform is an injective (one-to-one) map on the space of")
    print("Schwartz functions. The only function with a zero Fourier transform is the")
    print("zero function itself. This can be seen from the inverse Fourier transform:")
    print("  f(x) = Integral(\u00e2\u0178(xi) * exp(2*\u03c0*I*x*xi), (xi, -oo, oo))")
    print("  f(x) = Integral(0 * exp(2*\u03c0*I*x*xi), (xi, -oo, oo)) = 0\n")

    print("Therefore, it follows that f(x) must be the zero function.")
    print("\nThe final equation is:")
    final_equation = Eq(Function('f')(x), 0)
    print(f"  {final_equation}")

    print("\nThe number in the final equation is:")
    print(f"  {final_equation.rhs}")


if __name__ == '__main__':
    prove_moment_problem_solution()