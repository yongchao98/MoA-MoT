import sympy

def demonstrate_moment_problem_proof():
    """
    Demonstrates symbolically that if all moments of a Schwartz function are zero,
    the function itself must be zero.
    """

    # --- Setup: Define symbols and the core idea ---
    print("--- The Moment Problem for Schwartz Functions ---")
    print("Problem: If f(x) is a Schwartz function and Integral(x^k * f(x) dx) = 0 for all k=0, 1, 2,..., does f(x) = 0?")
    print("\nPlan: We will analyze this using the Fourier Transform, F(xi), of f(x).")

    # Define symbolic variables for demonstration
    x, xi = sympy.symbols('x xi', real=True)
    k = sympy.Symbol('k', integer=True, nonnegative=True)
    I = sympy.I # The imaginary unit

    # --- Step 1: Relate Fourier Transform Derivatives to Moments ---
    print("\nThe Fourier Transform is F(xi) = Integral(f(x) * exp(-2*pi*i*x*xi) dx).")
    print("Its k-th derivative at xi=0, F^(k)(0), relates to the k-th moment of f(x).")
    print("F^(k)(0) = Integral( f(x) * [d^k/dxi^k(exp(-2*pi*i*x*xi))]_(xi=0) dx )")

    # --- Step 2: Differentiate the Fourier Kernel ---
    print("\nLet's analyze the derivative of the kernel term, exp(-2*pi*i*x*xi), with respect to xi.")
    kernel = sympy.exp(-2 * sympy.pi * I * x * xi)
    
    # Calculate the general k-th derivative and evaluate it at xi = 0
    # The derivative of exp(a*xi) w.r.t. xi is a*exp(a*xi).
    # The k-th derivative is a^k*exp(a*xi), where a = -2*pi*i*x.
    kernel_k_derivative = (-2 * sympy.pi * I * x)**k * kernel
    kernel_k_derivative_at_0 = kernel_k_derivative.subs(xi, 0)
    
    print(f"The k-th derivative of the kernel with respect to xi is: ({kernel_k_derivative})")
    print(f"Evaluating this derivative at xi = 0 gives: {kernel_k_derivative_at_0}")
    
    # --- Step 3: Form the Key Equation ---
    print("\nSubstituting this back gives the full relationship:")
    # Using sympy's pretty print to format the equation
    moment = sympy.Integral(x**k * f(x), (x, -sympy.oo, sympy.oo))
    F_k_0 = sympy.Derivative(sympy.Function('F')(xi), (xi, k)).subs(xi, 0)
    
    equation = sympy.Eq(F_k_0, (-2 * sympy.pi * I)**k * moment)
    print("The k-th derivative of the Fourier transform at the origin is proportional to the k-th moment:")
    sympy.pprint(equation, use_unicode=True)
    
    # --- Step 4: Apply the Condition and Conclude ---
    print("\n--- Conclusion ---")
    print("We are given that the moment, Integral(x^k * f(x) dx), is 0 for all k.")
    print("Therefore, from the equation above, the k-th derivative of F(xi) at xi=0 must also be 0 for all k.")
    
    print("\nThe Taylor series for F(xi) around xi=0 is Sum[ (F^(k)(0) / k!) * xi^k ].")
    print("Since all derivatives F^(k)(0) are 0, the Taylor series is 0.")
    print("For a Schwartz function, its Fourier transform is analytic, so F(xi) must be identically 0 for all xi.")
    print("\nFinally, by the uniqueness of the inverse Fourier transform, if F(xi) = 0, then f(x) must be 0.")

demonstrate_moment_problem_proof()