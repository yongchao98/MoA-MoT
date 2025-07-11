import sympy

# Let's define the symbols to represent the physical quantities.
# We use sympy to display the mathematical expressions nicely.
gamma = sympy.Symbol('γ')  # Surface tension
r = sympy.Symbol('r')      # Radial position
xi = sympy.Function('ξ')(r) # Interfacial displacement as a function of r

def solve_and_display():
    """
    This function derives and displays the coefficients A(r) and B(r) for the governing equation.
    The derivation is based on the pressure balance between surface tension and electrostatic forces
    at a fluid interface in a cylindrical coordinate system, under a linear approximation.
    """
    
    # The governing equation is given in the form:
    # A(r) * d²ξ/dr² + B(r) * dξ/dr + C(r, ξ(r)) = 0
    
    # From the physical derivation (linearized Young-Laplace equation), the surface tension
    # part of the equation is γ * (d²ξ/dr² + (1/r) * dξ/dr).
    
    # By expanding this term: γ*d²ξ/dr² + (γ/r)*dξ/dr, we can identify
    # the coefficients of the second and first derivative terms.
    
    # Coefficient of the second derivative term, d²ξ/dr²
    A_r = gamma
    
    # Coefficient of the first derivative term, dξ/dr
    B_r = gamma / r
    
    # The term C(r, ξ(r)) represents the contribution from the electrostatic pressure,
    # which is not asked for.
    
    print("The governing linear differential equation for the interfacial shape ξ(r) is of the form:")
    print("A(r) * d²ξ/dr² + B(r) * dξ/dr + C(r, ξ(r)) = 0\n")
    
    print("Based on the derivation from the balance of forces, the coefficients are:")
    
    # Output the expression for A(r)
    print("A(r) =", sympy.pretty(A_r))
    
    # Output the expression for B(r)
    print("B(r) =", sympy.pretty(B_r))

    print("\nWhere:")
    print("γ (gamma) is the surface tension between the two fluids.")
    print("r is the radial coordinate.")

if __name__ == "__main__":
    solve_and_display()
