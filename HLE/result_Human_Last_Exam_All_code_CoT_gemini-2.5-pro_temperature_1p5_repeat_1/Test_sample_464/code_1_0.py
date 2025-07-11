import sympy

def solve_schwartz_moment_problem():
    """
    Shows that a Schwartz function with all zero moments must be the zero function.
    This is a logical demonstration using sympy for symbolic representation.
    """

    # Define the symbolic variables and functions
    x, xi = sympy.symbols('x xi', real=True)
    k = sympy.Symbol('k', integer=True, nonnegative=True)
    f = sympy.Function('f')
    F = sympy.Function('F') # Represents the Fourier transform of f

    print("--- Step 1: State the Premise ---")
    print("Let f(x) be a Schwartz class function.")
    print("The given condition is that all its moments are zero for k = 0, 1, 2, ...:")
    moment_condition = sympy.Eq(sympy.Integral(x**k * f(x), (x, -sympy.oo, sympy.oo)), 0)
    print(moment_condition)
    print("\n" + "="*50 + "\n")

    print("--- Step 2: Introduce the Fourier Transform ---")
    print("Let F(xi) be the Fourier transform of f(x).")
    fourier_def = sympy.Eq(F(xi), sympy.Integral(f(x) * sympy.exp(-2 * sympy.pi * sympy.I * x * xi), (x, -sympy.oo, sympy.oo)))
    print(fourier_def)
    print("\n" + "="*50 + "\n")

    print("--- Step 3: Relate Moments to Derivatives of the Transform ---")
    print("The k-th derivative of F(xi) can be found by differentiating under the integral sign:")
    derivative_of_F = sympy.Derivative(F(xi), (xi, k))
    moment_relation_rhs = sympy.Integral(f(x) * sympy.Derivative(sympy.exp(-2 * sympy.pi * sympy.I * x * xi), (xi, k)), (x, -sympy.oo, sympy.oo))
    moment_relation = sympy.Eq(derivative_of_F, moment_relation_rhs)
    print("Initially:")
    print(moment_relation)

    evaluated_derivative = sympy.Eq(derivative_of_F, sympy.Integral(f(x) * (-2 * sympy.pi * sympy.I * x)**k * sympy.exp(-2 * sympy.pi * sympy.I * x * xi), (x, -sympy.oo, sympy.oo)))
    print("\nAfter differentiation:")
    print(evaluated_derivative)
    print("\n" + "="*50 + "\n")


    print("--- Step 4: Evaluate at xi = 0 ---")
    print("Evaluating the k-th derivative at xi = 0 (since exp(0) = 1):")
    derivative_at_0 = derivative_of_F.subs(xi, 0)
    moment_relation_at_0_rhs = sympy.Integral((-2 * sympy.pi * sympy.I * x)**k * f(x), (x, -sympy.oo, sympy.oo))
    moment_relation_at_0 = sympy.Eq(derivative_at_0, moment_relation_at_0_rhs)
    print(moment_relation_at_0)

    print("\nFactoring out the constant (-2*pi*i)^k from the integral:")
    factored_moment = sympy.Eq(derivative_at_0, (-2 * sympy.pi * sympy.I)**k * sympy.Integral(x**k * f(x), (x, -sympy.oo, sympy.oo)))
    print(factored_moment)
    print("\n" + "="*50 + "\n")

    print("--- Step 5: Apply the Given Condition ---")
    print("We are given that the moment integral is 0. Substituting this into our equation:")
    # The equation has Integral(...) on the right side, so we substitute it with 0
    result_for_derivatives = sympy.Eq(derivative_at_0, (-2 * sympy.pi * sympy.I)**k * 0)
    print(result_for_derivatives)
    
    final_derivative_result = sympy.Eq(derivative_at_0, 0)
    print("\nThis simplifies to:")
    print(final_derivative_result)
    print("\nThis holds for all k = 0, 1, 2, ... So, all derivatives of F(xi) at xi=0 are zero.")
    print("\n" + "="*50 + "\n")

    print("--- Step 6: The Final Conclusion ---")
    print("A key theorem states that if a Schwartz function (like F(xi)) has all derivatives")
    print("at a single point equal to zero, the function must be identically zero.")
    conclusion_F = sympy.Eq(F(xi), 0)
    print(f"\nTherefore, we must have: {conclusion_F}")
    
    print("\nFinally, the Fourier transform is injective on Schwartz space. This means")
    print("if the transform F(xi) is the zero function, the original function f(x)")
    print("must also be the zero function.")

    final_conclusion = sympy.Eq(f(x), 0)
    print("\nSo, the final equation is:")
    print(final_conclusion)

solve_schwartz_moment_problem()