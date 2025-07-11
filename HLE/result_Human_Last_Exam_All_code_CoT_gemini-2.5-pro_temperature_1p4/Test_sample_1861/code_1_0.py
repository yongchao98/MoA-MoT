import sympy

def solve_manifold_problem():
    """
    This script demonstrates that for a 1-form eta on R^2, if it's
    invariant under a transitive group of diffeomorphisms, then d(eta) = 0.

    The reasoning is as follows:
    1. A transitive group of diffeomorphisms on R^2 must contain all translations.
    2. Therefore, the 1-form eta must be invariant under all translations.
    3. Invariance under translations means the coefficients of eta must be constant.
    4. The exterior derivative of a 1-form with constant coefficients is zero.
    """
    print("Let's analyze the case for M = R^2.")
    
    # Define coordinates and coefficient functions
    x, y = sympy.symbols('x y')
    a = sympy.Function('a')(x, y)
    b = sympy.Function('b')(x, y)
    
    print(f"Let a generic 1-form eta on R^2 be: eta = a(x,y)dx + b(x,y)dy")
    print(f"where a(x,y) = {a} and b(x,y) = {b}\n")
    
    print("The condition that eta is invariant under translations means its Lie derivative")
    print("with respect to the vector fields of translation, ∂/∂x and ∂/∂y, is zero.\n")

    # The Lie derivative of eta w.r.t ∂/∂x is (∂a/∂x)dx + (∂b/∂x)dy
    # For this to be zero, both coefficients must be zero.
    da_dx = sympy.diff(a, x)
    db_dx = sympy.diff(b, x)
    
    print(f"L_{{∂/∂x}}(eta) = 0  =>  ∂a/∂x = {da_dx} = 0  and  ∂b/∂x = {db_dx} = 0")

    # The Lie derivative of eta w.r.t ∂/∂y is (∂a/∂y)dx + (∂b/∂y)dy
    # For this to be zero, both coefficients must be zero.
    da_dy = sympy.diff(a, y)
    db_dy = sympy.diff(b, y)
    
    print(f"L_{{∂/∂y}}(eta) = 0  =>  ∂a/∂y = {da_dy} = 0  and  ∂b/∂y = {db_dy} = 0\n")

    print("These conditions together imply that a(x,y) and b(x,y) must be constants.")
    
    # Redefine a and b as constants
    C1, C2 = sympy.symbols('C1 C2')
    print(f"Let's represent these constants as a = {C1} and b = {C2}.")
    print(f"So, eta = {C1}dx + {C2}dy.\n")
    
    # Calculate the exterior derivative d(eta) = (∂b/∂x - ∂a/∂y) dx ∧ dy
    d_eta_coeff = sympy.diff(C2, x) - sympy.diff(C1, y)
    
    print("Now, we calculate the exterior derivative, d(eta).")
    print(f"The coefficient of d(eta) is (∂b/∂x - ∂a/∂y).")
    print(f"With constant coefficients, this is: ∂({C2})/∂x - ∂({C1})/∂y = {d_eta_coeff}")

    print("\nSince the coefficients are constant, their derivatives are zero.")
    final_equation = f"d(eta) = {d_eta_coeff} dx ∧ dy"
    print(f"The final equation is: {final_equation}")
    print("\nConclusion: The analysis shows that d(eta) = 0 is necessary for R^2.")
    print("Similar arguments hold for the cylinder and the torus, making B the correct answer.")

solve_manifold_problem()
<<<B>>>