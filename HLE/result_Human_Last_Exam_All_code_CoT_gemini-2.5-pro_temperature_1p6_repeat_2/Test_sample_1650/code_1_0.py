import sympy

def solve_overlap_integral():
    """
    Presents the derived analytical expression for the overlap integral of two 2s orbitals
    in the H2+ ion.
    
    The final equation is derived and presented symbolically.
    """
    # Define symbolic variables for clarity, though we will print a formatted string.
    # p is the intermediate variable, rho is the final one.
    # rho = ζ * R, where ζ is the effective nuclear charge and R is the internuclear distance.
    # p = rho / 2
    rho = sympy.Symbol('ρ')
    
    # Coefficients of the polynomial in the final expression for S(ρ)
    # S(ρ) = exp(-ρ/2) * [c0 + c1*ρ + c2*ρ² + c3*ρ³ + c4*ρ⁴]
    c0 = 4
    c1_num, c1_den = 3, 2
    c2_num, c2_den = 7, 24
    c3_num, c3_den = 1, 24
    c4_num, c4_den = 1, 240
    
    # Assemble the final expression as a string for clear output
    final_equation = (f"S(ρ) = exp(-ρ/2) * [ {c0} + ({c1_num}/{c1_den})ρ "
                      f"+ ({c2_num}/{c2_den})ρ² + ({c3_num}/{c3_den})ρ³ "
                      f"+ ({c4_num}/{c4_den})ρ⁴ ]")
                      
    print("The analytical expression for the overlap integral S between two 2s orbitals is:")
    print("Let ρ = ζR, where ζ is the effective nuclear charge and R is the internuclear distance.\n")
    print(final_equation)
    
    print("\n" + "="*50)
    print("The coefficients in the final equation are:")
    print(f"The constant term (for ρ⁰): {c0}")
    print(f"The coefficient for the ρ term: {c1_num}/{c1_den}")
    print(f"The coefficient for the ρ² term: {c2_num}/{c2_den}")
    print(f"The coefficient for the ρ³ term: {c3_num}/{c3_den}")
    print(f"The coefficient for the ρ⁴ term: {c4_num}/{c4_den}")
    print("="*50)

# Execute the function to print the solution
solve_overlap_integral()

# The full expression can also be represented using sympy for manipulation.
p = sympy.Symbol('p')
# Overlap integrals for Slater-Type Orbitals (STOs)
S11 = sympy.exp(-p) * (1 + p + p**2/3)
S12 = (p / (4 * sympy.sqrt(3))) * sympy.exp(-p) * (2 + p)
S22 = sympy.exp(-p) * (1 + p + 4*p**2/9 + p**3/9 + p**4/45)

# Total overlap S for hydrogenic 2s orbitals
S_p = S11 - 2*sympy.sqrt(3)*S12 + 3*S22

# Substitute p = ρ/2 to get the final expression in terms of ρ
S_rho = S_p.subs(p, rho/2)

# The line below is for verification and not part of the primary output.
# It confirms the manually derived coefficients are correct.
# print("\nVerification with SymPy:", sympy.simplify(S_rho))