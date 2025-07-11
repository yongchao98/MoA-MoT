import sympy
from sympy import I, pi, oo, exp, pretty

def solve_task():
    """
    This function demonstrates that if a Schwartz function has all zero moments,
    it must be the zero function. It does so by walking through the
    mathematical proof symbolically using the sympy library.
    """
    # Setup symbolic variables for the demonstration
    x, xi = sympy.symbols('x xi', real=True)
    k = sympy.Symbol('k', integer=True, nonnegative=True)
    f = sympy.Function('f')(x)

    print("--- The Moment Problem for Schwartz Functions ---")
    print("\nThe question is: If f is a Schwartz function such that all its moments are zero, is f necessarily the zero function?")
    print("Let's walk through the proof step-by-step.\n")

    print("Step 1: The Setup")
    print("Let f(x) be a Schwartz class function. A Schwartz function is infinitely differentiable")
    print("and decays faster than any polynomial at infinity, along with all its derivatives.")
    print("The k-th moment of f(x) is defined as M_k = integral(x**k * f(x) dx) over the real line.")
    print("The given condition is that M_k = 0 for all k = 0, 1, 2, ...")

    print("\nStep 2: The Fourier Transform and its Derivatives")
    print("The Fourier transform of f(x) is defined as F(xi) = integral(f(x) * exp(-2*pi*I*x*xi) dx).")
    print("Let's find the relationship between the derivatives of F(xi) at xi=0 and the moments M_k.")
    
    # The exponential term in the Fourier transform integral
    integrand_exp_term = exp(-2*I*pi*x*xi)
    
    # Let's compute the k-th derivative of this term with respect to xi
    deriv_exp_term = sympy.diff(integrand_exp_term, xi, k)
    
    print(f"\nThe k-th derivative of exp(-2*pi*I*x*xi) with respect to xi is:")
    print(pretty(deriv_exp_term, use_unicode=False))
    
    # Now, let's evaluate this derivative at xi = 0
    deriv_exp_term_at_0 = deriv_exp_term.subs(xi, 0)
    
    print(f"\nEvaluating this derivative at xi=0 gives:")
    print(pretty(deriv_exp_term_at_0, use_unicode=False))
    
    print("\nBecause f(x) is a Schwartz function, we can legally swap differentiation and integration.")
    print("So, the k-th derivative of F(xi) evaluated at xi=0, denoted F^(k)(0), is:")
    print("F^(k)(0) = integral(f(x) * (" + str(pretty(deriv_exp_term_at_0, use_unicode=False)) + ") dx)")
    
    # The constant factor that can be pulled out of the integral
    constant_factor_symbolic = (-2*I*pi)**k
    
    print("\nThis gives us the crucial relationship. First, we factor out the constant terms:")
    # Print the equation showing the relationship. This satisfies the "output each number" requirement
    # in a symbolic context by printing the constituent parts of the equation.
    equation_part1 = f"F^(k)(0) = {pretty(constant_factor_symbolic, use_unicode=False)} * integral(x**k * f(x) dx)"
    print(equation_part1)
    
    print("\nRecognizing the integral as the k-th moment M_k, the equation becomes:")
    equation_part2 = f"F^(k)(0) = {pretty(constant_factor_symbolic, use_unicode=False)} * M_k"
    print(equation_part2)

    print("\nStep 3: Applying the Condition")
    print("We are given that M_k = 0 for all non-negative integers k.")
    print(f"Substituting M_k = 0 into our equation: F^(k)(0) = {pretty(constant_factor_symbolic, use_unicode=False)} * 0 = 0.")
    print("This proves that all derivatives of the Fourier transform at xi=0 are zero.")

    print("\nStep 4: The Conclusion using Analyticity")
    print("A key property of Schwartz functions is that their Fourier transforms are analytic functions.")
    print("An analytic function whose value and all derivatives at a single point are zero must be the zero function.")
    print("Since F^(k)(0) = 0 for all k, the Taylor series of F(xi) around xi=0 is identically zero.")
    print("This implies that F(xi) = 0 for all xi.")
    
    print("\nFinally, the Fourier transform on the space of Schwartz functions is an invertible mapping.")
    print("Therefore, if F(xi) is the zero function, its inverse Fourier transform f(x) must also be the zero function.")
    
    print("\n-------------------")
    print("Final Answer: Yes, it follows that f = 0.")
    print("-------------------")

solve_task()