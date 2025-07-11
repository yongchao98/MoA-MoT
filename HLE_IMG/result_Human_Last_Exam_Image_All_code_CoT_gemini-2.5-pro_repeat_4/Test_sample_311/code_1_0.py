import sympy
import math
import numpy as np

def calculate_Ci():
    """
    Calculates the list of C_i coefficients based on properties derived from each visualization.
    This function encapsulates the logic from steps 1-3 of the plan.
    """
    # Group is D6. Order=12, Irreps=6, Exponent=6. Char table sum=8.
    
    # R_ji values for each V_i. {R1, R2, R3, R4}
    # R1: Char table sum, R2: #irreps, R3: order, R4: exponent
    R = {
        1: [0, 0, 12, 6],  # V1 (Cayley Graph)
        2: [0, 0, 12, 6],  # V2 (Cycle Graph)
        3: [0, 0, 0, 0],   # V3 (Anomalous graph, no info on D6)
        4: [0, 0, 12, 0],  # V4 (Resistance Matrix)
        5: [8, 6, 12, 6],  # V5 (Character Table)
        6: [0, 6, 12, 6],  # V6 (Multiplication Table)
        7: [0, 0, 12, 6],  # V7 (Adj Matrix of Cycle Graph)
        8: [0, 0, 12, 6],  # V8 (Incidence Matrix of Cycle Graph)
        9: [0, 0, 12, 6],  # V9 (Kirchhoff Matrix of Cycle Graph)
    }
    
    C = []
    for i in range(1, 10):
        vals = R[i]
        sum_vals = sum(vals)
        sum_sq_vals = sum(v**2 for v in vals)
        
        if sum_vals == 0:
            chm = 0
        else:
            chm = sum_sq_vals / sum_vals
            
        C.append(math.floor(chm))
        
    return C

def solve_problem():
    """
    Solves the entire problem as laid out in the plan.
    """
    # Step 3: Compute C_i coefficients
    C = calculate_Ci()
    # C = [10, 10, 0, 12, 8, 9, 10, 10, 10]

    # Step 4: Construct Polynomials Q(x) and S(x)
    x = sympy.Symbol('x')
    
    # Coefficients for Q(x) = -C2*x^2 + C4*x^4 - C6*x^6 + C8*x^8
    Q_coeffs = {-2: -C[1], 4: C[3], -6: -C[5], 8: C[7]}
    Q = -C[1]*x**2 + C[3]*x**4 - C[5]*x**6 + C[7]*x**8
    
    # Coefficients for S(x) = C1*x - C3*x^3 + C5*x^5 - C7*x^7 + C9*x^9
    S_coeffs = {1: C[0], -3: -C[2], 5: C[4], -7: -C[6], 9: C[8]}
    S = C[0]*x - C[2]*x**3 + C[4]*x**5 - C[6]*x**7 + C[8]*x**9

    # Step 5: Calculate Tr(M1)
    M1_matrix = sympy.polys.matrices.bezout(Q, S, x)
    trace_M1 = M1_matrix.to_dense().trace()

    # Step 5: Calculate Tr(M2)
    P1 = Q
    P2 = x**10 + S
    
    p_coeffs = P1.as_poly(x).all_coeffs()
    q_coeffs = P2.as_poly(x).all_coeffs()
    
    m = sympy.degree(P1, x)
    n = sympy.degree(P2, x)
    
    p_m = p_coeffs[0] if p_coeffs else 0
    q_0 = q_coeffs[-1] if q_coeffs and P2.as_poly(x).TC() == 0 else 0
    
    trace_M2 = n * p_m + m * q_0

    # Step 6: Final Calculation
    # T = Tr(M1 tensor I2 + M2) = Tr(M1)*Tr(I2) + Tr(M2)
    trace_I2 = 2
    T = trace_M1 * trace_I2 + trace_M2

    print("The final equation is T = 2 * Tr(M1) + Tr(M2)")
    print(f"Tr(M1) = {trace_M1}")
    print(f"Tr(M2) = {trace_M2}")
    print(f"T = 2 * ({trace_M1}) + {trace_M2} = {T}")
    print(f"<<<{T}>>>")

solve_problem()