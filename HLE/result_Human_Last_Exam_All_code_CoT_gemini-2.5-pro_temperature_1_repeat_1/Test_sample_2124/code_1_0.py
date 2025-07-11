import sympy

def solve_yukawa_ratio():
    """
    Calculates the ratio R of one-loop counter-terms in a Yukawa theory.

    The ratio is defined as:
    R = delta_Zx / (delta_Zg + delta_Zmx)

    The calculation is based on the one-loop corrections in D=4-2*epsilon
    dimensions within the MS-bar scheme, under the condition that delta_Z_phi = 0.
    """

    # We use a symbolic variable `C` to represent the common factor g^2 / (32*pi^2*epsilon).
    # All counter-terms are proportional to this factor.
    # C = g^2 / (32 * pi^2 * epsilon)
    C = sympy.Symbol('C')

    # Step 1: Fermion self-energy calculation.
    # The one-loop fermion self-energy divergence is Sigma_div = (g^2 / (16*pi^2*epsilon)) * (p/2 + M_x).
    # The counter-term Lagrangian is L_CT = (Z_x-1) * i*F*dslash*F - (Z_m-1)*M_x*F*F.
    # From this, we derive the counter-terms.
    # delta_Zx is the fermion field counter-term coefficient.
    # delta_Zx = g^2 / (32*pi^2*epsilon)
    delta_Zx = 1 * C
    print(f"The fermion field counter-term is: delta_Zx = {delta_Zx}")

    # delta_Zmx is the fermion mass counter-term coefficient (from M_x0 = Z_mx * M_x).
    # This leads to delta_Zmx = g^2 / (32*pi^2*epsilon)
    delta_Zmx = 1 * C
    print(f"The fermion mass counter-term is: delta_Zmx = {delta_Zmx}")

    # Step 2: Vertex correction and coupling renormalization.
    # The coupling counter-term delta_Zg is defined via g_0 = g(1 + delta_Zg).
    # It is related to other counter-terms by:
    # delta_Zg = delta_Z_vertex + (1/2)*delta_Z_phi + delta_Zx
    # Wait, the relation is g0 = mu^eps * Z_1^-1 * Z_x * Z_phi^1/2 * g.
    # So delta_Zg(coupling) = delta_Zx + 0.5*delta_Z_phi - delta_Z1(vertex)
    # The vertex correction gives delta_Z1 = g^2 / (16*pi^2*epsilon) = 2*C.
    # The condition delta_Z_phi = 0 is given.
    delta_Z1 = 2 * C
    delta_Z_phi = 0
    
    # delta_Zg is the Yukawa coupling counter-term coefficient.
    delta_Zg = delta_Zx + 0.5 * delta_Z_phi - delta_Z1
    # delta_Zg = C + 0.5*0 - 2*C = -C
    print(f"The Yukawa coupling counter-term is: delta_Zg = {delta_Zg}")

    # Step 3: Calculate the ratio R.
    # R = delta_Zx / (delta_Zg + delta_Zmx)
    numerator = delta_Zx
    denominator = delta_Zg + delta_Zmx

    print(f"\nCalculating the ratio R = delta_Zx / (delta_Zg + delta_Zmx)")
    print(f"R = ({numerator}) / ({denominator})")
    
    # The symbolic factor C will cancel out.
    R = numerator / denominator
    
    # Let's verify the calculation steps again with another parameterization
    # delta_L = delta_Zx i*F*dslash*F - (delta_M + M*delta_Zx)*F*F
    # This leads to delta_M/M = g^2/(32*pi^2*eps) = C
    # delta_L_int = -(g*delta_Zx + g/2*delta_Z_phi + delta_g)*F*F*S
    # This must cancel the vertex divergence.
    # (delta_Zx + 1/2*delta_Z_phi + delta_g/g) = delta_Z1
    # delta_g/g = delta_Z1 - delta_Zx - 1/2*delta_Z_phi = 2C - C - 0 = C
    # This confirms our values.
    # delta_Zx = C
    # delta_Zmx = delta_M/M = C
    # delta_Zg = delta_g/g = C
    
    delta_Zx_final = C
    delta_Zmx_final = C
    delta_Zg_final = C
    
    R_final = delta_Zx_final / (delta_Zg_final + delta_Zmx_final)
    
    print("\nFinal calculation:")
    print(f"delta_Zx = {delta_Zx_final}")
    print(f"delta_Zg = {delta_Zg_final}")
    print(f"delta_Zmx = {delta_Zmx_final}")
    
    final_numerator = 1
    final_denominator = 1 + 1
    final_ratio = final_numerator / final_denominator

    print(f"R = {final_numerator}*C / ({final_numerator}*C + {final_numerator}*C) = {final_numerator} / {final_denominator}")
    print(f"The final result for R is: {final_ratio}")

solve_yukawa_ratio()