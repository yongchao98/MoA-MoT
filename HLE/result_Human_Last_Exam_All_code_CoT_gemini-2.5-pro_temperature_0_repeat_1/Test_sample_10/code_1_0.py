import sympy

def solve():
    """
    Analyzes the derivatives of four activation functions to find the one
    whose derivative cannot be written as a simple polynomial of the sigmoid function.
    """
    # Define symbolic variables
    x, beta = sympy.symbols('x beta')

    # --- Introduction ---
    print("Analyzing the first derivative of each function to check its relationship with the sigmoid function sigma(u) = 1/(1 + exp(-u)).")
    print("The key property of the sigmoid derivative is sigma'(u) = sigma(u) * (1 - sigma(u)).")
    print("We will check if each T'(x) can be written as a polynomial of a sigmoid function.")
    print("="*70)

    # --- Analysis of T3 (Softplus) ---
    print("Analysis of T3(x) = log(1 + exp(x))")
    T3 = sympy.log(1 + sympy.exp(x))
    T3_prime = sympy.diff(T3, x)
    # T3_prime simplifies to exp(x)/(1 + exp(x)), which is exactly sigma(x)
    s_x = sympy.Symbol('sigma(x)')
    print(f"The derivative is T3'(x) = {sympy.simplify(T3_prime)}")
    print(f"This is exactly the sigmoid function. It can be written as a simple polynomial of degree 1 in sigma(x):")
    print(f"T3'(x) = {s_x}")
    print("Result: T3' has a direct connection to sigma(x).")
    print("="*70)

    # --- Analysis of T1 (Swish) ---
    print("Analysis of T1(x) = x / (1 + exp(-beta*x))")
    # T1 is x * sigma(beta*x)
    s_beta = sympy.Symbol('sigma(beta*x)')
    # Using the product rule and the sigmoid derivative property:
    # T1'(x) = d/dx(x * sigma(beta*x)) = sigma(beta*x) + x * sigma'(beta*x) * beta
    # T1'(x) = sigma(beta*x) + beta*x * sigma(beta*x) * (1 - sigma(beta*x))
    T1_prime_in_sigma = s_beta + beta * x * s_beta * (1 - s_beta)
    print("The derivative can be written as a polynomial in sigma(beta*x):")
    print(f"T1'(x) = {s_beta} + {beta}*x*{s_beta}*(1 - {s_beta})")
    print(f"Expanding this gives: {sympy.expand(T1_prime_in_sigma)}")
    print("Result: T1' has a direct polynomial connection to sigma(beta*x).")
    print("="*70)

    # --- Analysis of T4 (GELU Approximation) ---
    print("Analysis of T4(x) = 0.5*x*(1 + tanh(sqrt(2/pi)*(x + 0.044715*x**3)))")
    # The term (1 + tanh(u)) can be rewritten as 2*sigma(2u).
    # So, T4(x) = 0.5*x*(2*sigma(2z)) = x*sigma(2z), where z = sqrt(2/pi)*(x + 0.044715*x**3)
    z = sympy.sqrt(2/sympy.pi) * (x + 0.044715 * x**3)
    z_prime = sympy.diff(z, x)
    s_2z = sympy.Symbol('sigma(2*z(x))')
    # Using the product rule similar to T1:
    # T4'(x) = sigma(2z) + x * d/dx(sigma(2z)) = sigma(2z) + x * sigma'(2z) * 2 * z'
    # T4'(x) = sigma(2z) + 2*x*z' * sigma(2z) * (1 - sigma(2z))
    T4_prime_in_sigma = s_2z + (2 * x * z_prime) * s_2z * (1 - s_2z)
    print("The derivative can be written as a polynomial in sigma(2*z(x)):")
    # We print the structure of the equation.
    print(f"T4'(x) = {s_2z} + (2*x*({sympy.simplify(z_prime)})) * {s_2z} * (1 - {s_2z})")
    print("Result: T4' has a direct polynomial connection to sigma(2*z(x)).")
    print("="*70)

    # --- Analysis of T2 (Mish) ---
    print("Analysis of T2(x) = x * ((-1 + (1 + exp(x))**2) / (1 + (1 + exp(x))**2))")
    # This is the Mish function, T2(x) = x * tanh(softplus(x)).
    # We will express its derivative in terms of s = sigma(x).
    s = sympy.Symbol('sigma(x)')
    # The derivative T2'(x) can be shown to be:
    # T2'(x) = f(s) + x * g(s), where f(s) and g(s) are rational functions of s.
    numerator1 = 2*s - s**2
    denominator1 = s**2 - 2*s + 2
    
    numerator2 = 4*s*(1-s)**2
    denominator2 = (s**2 - 2*s + 2)**2
    
    print("The derivative can be written in terms of s = sigma(x) as follows:")
    print(f"T2'(x) = ({numerator1}) / ({denominator1}) + x * (({numerator2}) / ({denominator2}))")
    print("This is a rational function of sigma(x), but it is not a polynomial because of the terms in the denominator.")
    print("Result: T2's connection to sigma(x) is not a simple polynomial one.")
    print("="*70)

    # --- Conclusion ---
    print("\nConclusion:")
    print("The derivatives of T1, T3, and T4 are all polynomials in a corresponding sigmoid function.")
    print("The derivative of T2 is a more complex rational function of the sigmoid function, not a polynomial.")
    print("Therefore, T2 is the function whose derivative is the most distinct and cannot be written as a (simple polynomial) function of the sigmoid.")

solve()