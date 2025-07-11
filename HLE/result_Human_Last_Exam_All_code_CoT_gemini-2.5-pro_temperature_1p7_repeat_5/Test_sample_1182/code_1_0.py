import sympy as sp

def solve_curve_reduction():
    """
    Calculates the number of double points in the stable reduction of the given curve above p=2.
    """
    x = sp.Symbol('x')
    
    # Original polynomial f(x)
    f_orig = 8*x**5 + 4*x**4 + 4*x**3 + x**2 + 8*x
    
    # Step 1: Calculate the genus of the original curve C: y^2 = f(x)
    # The degree of the polynomial f(x) is 5.
    deg_f = 5
    # The genus g is given by floor((deg(f) - 1) / 2)
    g = (deg_f - 1) // 2
    
    print(f"Step 1: The original curve is y^2 = 8*x**5 + 4*x**4 + 4*x**3 + 1*x**2 + 8*x.")
    print(f"The degree of the polynomial is {deg_f}.")
    print(f"The genus of the original curve is g = floor(({deg_f} - 1) / 2) = {g}.")
    print("-" * 20)

    # Step 2: The reduction modulo 2 of the original equation is y^2 = x^2.
    # This is a non-reduced line, which is an unstable reduction. We need to find a better model.
    print("Step 2: The reduction of this equation modulo 2 is y^2 = x^2, which is highly singular.")
    print("We need to find a better model by a change of variables.")
    print("-" * 20)

    # Step 3: Find a better model via coordinate change x = 2^k * x', y = 2^m * y'
    # The valuations of coefficients are: v2(8)=3, v2(4)=2, v2(1)=0.
    # After analyzing the valuations for x = 2^k * x', we find k=3 is the optimal choice.
    # This leads to m=3.
    # The transformation is x = 8x', y = 8y'.
    k=3
    m=3
    x_prime = sp.Symbol("x'")
    
    # New polynomial f'(x')
    f_new = (2**(12) * x_prime**5 + 2**8 * x_prime**4 + 2**5 * x_prime**3 + 
             x_prime**2 + x_prime)
    
    print(f"Step 3: We apply the coordinate change x = {2**k}*x' and y = {2**m}*y'.")
    print(f"The new model is y'^2 = {2**(12)}*x'^5 + {2**8}*x'^4 + {2**5}*x'^3 + 1*x'^2 + 1*x'.")
    print("-" * 20)
    
    # Step 4: Analyze the new model's reduction modulo 2
    # The reduction of f'(x') modulo 2
    f_reduced = x_prime**2 + x_prime
    
    # Check if the reduced curve y'^2 = x'^2 + x' is smooth.
    # A curve y^2 = h(x) is smooth if h(x) has no repeated roots.
    # h(x') = x'^2 + x' = x'(x'+1). Roots are 0 and 1, which are distinct in F_2.
    # Thus, the reduced curve is smooth.
    print("Step 4: We reduce the new model modulo 2.")
    print(f"The reduced equation is y'^2 = x'^2 + x'.")
    print("This curve is smooth because the polynomial x'^2 + x' has distinct roots (0 and 1) in the field F_2.")
    print("-" * 20)

    # Step 5: Calculate the genus of the reduced curve
    # The reduced curve is y'^2 = x'^2 + x', a plane conic.
    deg_f_reduced = 2
    # It is a smooth curve of degree 2, so its genus is g_s = (d-1)(d-2)/2 = 0.
    g_s = (deg_f_reduced - 1) * (deg_f_reduced - 2) // 2
    print(f"Step 5: Calculate the genus of the reduced curve y'^2 = x'^2 + x'.")
    print(f"This is a smooth curve of degree {deg_f_reduced}. Its genus g_s is {g_s}.")
    print("-" * 20)

    # Step 6: Calculate the number of double points
    # The number of double points delta is given by g - g_s
    delta = g - g_s
    
    print("Step 6: The number of double points (nodes) in the stable reduction is the difference between the genera.")
    print(f"Number of double points, δ = g - g_s")
    print(f"δ = {g} - {g_s} = {delta}")

    return delta

# Run the function and print the final answer
final_answer = solve_curve_reduction()
print(f"\nThe number of double points is {final_answer}.")
print(f"<<<{final_answer}>>>")
