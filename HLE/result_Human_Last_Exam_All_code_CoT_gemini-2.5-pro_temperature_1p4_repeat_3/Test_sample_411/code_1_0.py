import sympy as sp

def analyze_stbc():
    """
    Analyzes the diversity order of three space-time block codes (STBCs)
    and determines which provides the maximum diversity.
    """
    # Define symbolic variables for the differences of two complex symbols
    # dx1 = x1' - x1", dx2 = x2' - x2"
    dx1 = sp.Symbol('Δx1')
    dx2 = sp.Symbol('Δx2')
    
    print("--- Diversity Analysis of Space-Time Codes ---")
    print("The diversity order is the minimum rank of the codeword difference matrix ΔS.")
    print("We check this by calculating the determinant of ΔS.\n")

    # --- Code Sa ---
    print("(a) Analysis of Code Sa:")
    # Codeword difference matrix ΔSa
    delta_S_a = sp.Matrix([
        [dx1, dx2],
        [dx2, dx1]
    ])
    # Determinant of ΔSa
    det_S_a = sp.det(delta_S_a)
    print(f"For Sa = [[x1, x2], [x2, x1]], the determinant of the difference matrix is:")
    # The output needs to show each number in the final equation.
    # The equation is det = (Δx1)^2 - (Δx2)^2. Let's print it formatted.
    print(f"det(ΔSa) = ({dx1})**2 - ({dx2})**2 = {det_S_a}")
    print("This determinant is zero if Δx1 = Δx2 or Δx1 = -Δx2.")
    print("Since this can happen for non-zero symbol differences, the matrix can be rank-deficient.")
    print("Diversity Order for Sa = 1.\n")

    # --- Code Sb ---
    print("(a) Analysis of Code Sb:")
    # Codeword difference matrix ΔSb
    # We need to handle the conjugate properly. We define dx1c as the conjugate of dx1.
    dx1c = sp.Symbol('Δx1*') 
    delta_S_b = sp.Matrix([
        [dx1, dx2],
        [dx2, dx1c]
    ])
    # Determinant of ΔSb
    det_S_b = sp.det(delta_S_b)
    print(f"For Sb = [[x1, x2], [x2, x1*]], the determinant of the difference matrix is:")
    print(f"det(ΔSb) = ({dx1})*({dx1c}) - ({dx2})**2 = |{dx1}|**2 - ({dx2})**2")
    print("This determinant can be zero. For example, if Δx2 is real and |Δx1| = |Δx2|.")
    print("Since this can happen for non-zero symbol differences, the matrix can be rank-deficient.")
    print("Diversity Order for Sb = 1.\n")

    # --- Code Sc ---
    print("(a) Analysis of Code Sc:")
    # Codeword difference matrix ΔSc
    dx2c = sp.Symbol('Δx2*')
    delta_S_c = sp.Matrix([
        [-dx1c, dx2],
        [-dx2c, -dx1]
    ])
    # Determinant of ΔSc
    det_S_c = sp.det(delta_S_c)
    print(f"For Sc = [[-x1*, x2], [-x2*, -x1]], the determinant of the difference matrix is:")
    print(f"det(ΔSc) = (-{dx1c})*(-{dx1}) - ({dx2})*(-{dx2c}) = |{dx1}|**2 + |{dx2}|**2")
    print("This determinant is the sum of two squared magnitudes.")
    print("It is zero only if both Δx1=0 and Δx2=0, which is not allowed for distinct codewords.")
    print("Therefore, the matrix is always full rank.")
    print("Diversity Order for Sc = 2.\n")
    
    # --- Part (b) Conclusion ---
    print("--- Final Conclusion ---")
    print("(a) The diversity orders are:")
    print("    - Code Sa: 1")
    print("    - Code Sb: 1")
    print("    - Code Sc: 2")
    print("\n(b) Which code provides the maximum directivity?")
    print("Assuming 'directivity' refers to robustness against fading, this corresponds to the maximum diversity order.")
    print("Code Sc has the highest diversity order (2), which is the full diversity for a 2x1 system.")
    print("Therefore, Code Sc provides the maximum directivity.")

if __name__ == '__main__':
    analyze_stbc()