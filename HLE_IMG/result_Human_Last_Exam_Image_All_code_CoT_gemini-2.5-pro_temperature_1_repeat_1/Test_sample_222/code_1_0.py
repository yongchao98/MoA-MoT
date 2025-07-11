import sympy

def solve_complex_integral():
    """
    Calculates the imaginary part of the sum of the given contour integrals.
    """
    z = sympy.Symbol('z')
    
    # Poles and their total winding numbers N_k = Ind_C1(z_k) + Ind_C2(z_k)
    poles_info = {
        sympy.S(3)/2: 2,
        -1: 2,
        -2: 2,
        -3: 1,
    }
    
    # The function f(z)
    f = (z / (z - sympy.S(3)/2)) * sympy.gamma(z)
    
    # Calculate residues
    residues = {pole: sympy.residue(f, z, pole) for pole in poles_info}
    
    # Calculate the sum R = Sum(N_k * Res(f, z_k))
    R_terms = {pole: poles_info[pole] * residues[pole] for pole in poles_info}
    R = sum(R_terms.values())
    
    # The imaginary part of the integral sum is 2 * pi * R
    imaginary_part = 2 * sympy.pi * R
    
    # Print the step-by-step calculation
    print("The imaginary part of the sum of the integrals is Im(I) = 2 * pi * R, where R is the sum of residues multiplied by their total winding numbers.")
    print("R = N(1.5)*Res(f, 1.5) + N(-1)*Res(f, -1) + N(-2)*Res(f, -2) + N(-3)*Res(f, -3)")
    print("\nWith calculated values:")
    
    res_1_5 = residues[sympy.S(3)/2]
    res_m1 = residues[-1]
    res_m2 = residues[-2]
    res_m3 = residues[-3]

    print(f"Res(f, 3/2) = {sympy.pretty(res_1_5)}")
    print(f"Res(f, -1) = {res_m1}")
    print(f"Res(f, -2) = {res_m2}")
    print(f"Res(f, -3) = {res_m3}\n")
    
    # Print the equation with all numbers
    print(f"Im(I) = 2 * pi * ( ({poles_info[sympy.S(3)/2]} * ({res_1_5})) + ({poles_info[-1]} * ({res_m1})) + ({poles_info[-2]} * ({res_m2})) + ({poles_info[-3]} * ({res_m3})) )")
    
    term_1_5 = R_terms[sympy.S(3)/2]
    term_m1 = R_terms[-1]
    term_m2 = R_terms[-2]
    term_m3 = R_terms[-3]
    
    print(f"Im(I) = 2 * pi * ( {term_1_5} + {term_m1} + {term_m2} + {term_m3} )")
    
    # Sum the rational parts
    rational_sum = term_m1 + term_m2 + term_m3
    
    print(f"Im(I) = 2 * pi * ( {rational_sum} + {term_1_5} )")
    
    # Final symbolic expression
    final_expression = 2 * sympy.pi * rational_sum + 2 * sympy.pi * term_1_5
    print(f"Im(I) = {final_expression}")
    
    # Final numerical value
    print(f"\nNumerical value: {imaginary_part.evalf()}")

solve_complex_integral()
<<<14.5590321208948>>>