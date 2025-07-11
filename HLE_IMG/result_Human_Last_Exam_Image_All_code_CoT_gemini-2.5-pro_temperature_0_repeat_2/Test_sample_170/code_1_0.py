import numpy as np
from scipy.optimize import fsolve

def solve_hvac_system():
    """
    Formulates and solves the non-linear optimization problem for the HVAC system.
    """
    # Given parameters in per-unit (p.u.), assuming S_base = 100 MVA
    R_S = 0.02  # p.u.
    X_S = 0.10  # p.u.
    Z_S_sq = R_S**2 + X_S**2
    V_A = 1.0  # p.u.
    V_B1 = 0.85  # p.u.
    V_B2 = 1.0  # p.u.
    S_base = 100.0  # MVA

    # The problem is defined by a system of 5 non-linear equations
    # for 5 variables: Q (Q_comp), c1 (cos(d1)), s1 (sin(d1)), c2 (cos(d2)), s2 (sin(d2))
    def equations(vars):
        Q, c1, s1, c2, s2 = vars
        
        # Equation 1: From real power balance, eliminating P_L and Q_L
        # 0.10 * Q = 0.85*c1 - c2 + 0.2775
        eq1 = X_S * Q - (V_A * V_B1 * c1 - V_B1**2) + (V_A * V_B2 * c2 - V_B2**2)
        
        # Equation 2: From reactive power balance, eliminating P_L and Q_L
        # 0.02 * Q = 0.85*s1 - s2
        eq2 = R_S * Q - (V_A * V_B1 * s1) + (V_A * V_B2 * s2)
        
        # Equation 3: From the 4% loss increase constraint
        # 2(1 - c2) = 1.04 * (1.7225 - 1.7*c1)
        term1 = V_A**2 + V_B2**2 - 2 * V_A * V_B2 * c2
        term2 = V_A**2 + V_B1**2 - 2 * V_A * V_B1 * c1
        eq3 = term1 - 1.04 * term2
        
        # Equation 4: Trigonometric identity for angle d1
        eq4 = c1**2 + s1**2 - 1
        
        # Equation 5: Trigonometric identity for angle d2
        eq5 = c2**2 + s2**2 - 1
        
        return [eq1, eq2, eq3, eq4, eq5]

    # Initial guess for the solver [Q, cos(d1), sin(d1), cos(d2), sin(d2)]
    # Assume Q is positive, angles are small and negative
    initial_guess = [0.5, 0.95, -0.3, 0.98, -0.2]
    
    # Solve the system of equations
    solution, infodict, ier, mesg = fsolve(equations, initial_guess, full_output=True)

    if ier == 1:
        Q_opt_pu = solution[0]
        c2 = solution[3]
        
        # Convert results to standard units
        Q_opt_mvar = Q_opt_pu * S_base
        
        # Calculate final system real power loss in MW
        P_loss_final_pu = (V_A**2 + V_B2**2 - 2 * V_A * V_B2 * c2) / Z_S_sq * R_S
        P_loss_final_mw = P_loss_final_pu * S_base
        
        print("Optimization Results:")
        print(f"The minimum reactive power injection required is Q_opt = {Q_opt_mvar:.2f} MVAR.")
        print(f"The final system real power loss is P_loss = {P_loss_final_mw:.2f} MW.")
        
        # For completeness, let's show the final equation for Q_opt
        print("\nThis was found by solving the system of equations, which can be reduced to an equation for Q_opt.")
        print("A key derived relationship is:")
        print(f"{X_S:.2f} * Q_opt = {V_A*V_B1:.2f}*cos(d1) - cos(d2) + {V_B2**2 - V_B1**2:.4f}")
        print(f"where d1 and d2 are constrained by the power flow and the 4% loss increase.")
        print(f"The solved value for Q_opt is {Q_opt_mvar:.2f} MVAR.")

    else:
        print("A solution could not be found with the given parameters.")
        print("This indicates that the problem statement may contain inconsistent values,")
        print("as the physical constraints cannot all be met simultaneously.")
        print(f"Solver message: {mesg}")

solve_hvac_system()
<<<42.35>>>