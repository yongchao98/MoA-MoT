def solve_critical_exponent():
    """
    Calculates the critical exponent p for the L^p integrability of the function I(a).
    """

    # Step 1: Define dimensions of the coefficient subspaces
    d_L = 2  # Dimension of the linear part (a1, a2)
    d_Q = 3  # Dimension of the quadratic part (a3, a4, a5)
    d_C = 4  # Dimension of the cubic part (a6, a7, a8, a9)
    n_vars = 2 # Number of integration variables (x, y)

    # Step 2: Define decay rates for each dominant region
    # The decay rate for an oscillatory integral with a homogeneous polynomial phase of degree d
    # in n variables is R^(-n/d).
    delta_L = 1.0  # For linear phase, decay is R^-1
    delta_Q = n_vars / 2.0  # For quadratic phase (d=2)
    delta_C = n_vars / 3.0  # For cubic phase (d=3)

    print("Analysis of the integral convergence:")
    print("-" * 35)

    # Step 3: Calculate the critical exponent for each region

    # Case 1: Cubic terms dominate (Region D_C)
    # The integral behaves like integral(R^(d_L+d_Q) * (R^(-delta_C))^p * R^(d_C-1) dR)
    # which is integral(R^(d_L+d_Q+d_C-1 - p*delta_C) dR)
    power_C = d_L + d_Q + d_C - 1
    # For convergence, power_C - p*delta_C < -1
    # p*delta_C > power_C + 1
    # p > (power_C + 1) / delta_C
    p_c_C = (power_C + 1) / delta_C
    print(f"In the region where cubic terms dominate (D_C):")
    print(f"  Dimension of linear space d_L = {d_L}")
    print(f"  Dimension of quadratic space d_Q = {d_Q}")
    print(f"  Dimension of cubic space d_C = {d_C}")
    print(f"  Decay rate delta_C = n_vars / 3 = {n_vars}/3 = {delta_C:.4f}")
    print(f"  The radial part of the integral behaves like R^({d_L} + {d_Q} + {d_C} - 1 - p * {delta_C:.2f}) = R^({power_C} - p*{delta_C:.2f})")
    print(f"  For convergence, {power_C} - p*{delta_C:.2f} < -1")
    print(f"  This implies p > ({power_C} + 1) / {delta_C:.4f} = {p_c_C}")
    print("-" * 35)

    # Case 2: Quadratic terms dominate (Region D_Q)
    # The integral behaves like integral(R^(d_L+d_C) * (R^(-delta_Q))^p * R^(d_Q-1) dR)
    power_Q = d_L + d_C + d_Q - 1
    # For convergence, power_Q - p*delta_Q < -1
    # p > (power_Q + 1) / delta_Q
    p_c_Q = (power_Q + 1) / delta_Q
    print(f"In the region where quadratic terms dominate (D_Q):")
    print(f"  Dimension of linear space d_L = {d_L}")
    print(f"  Dimension of cubic space d_C = {d_C}")
    print(f"  Dimension of quadratic space d_Q = {d_Q}")
    print(f"  Decay rate delta_Q = n_vars / 2 = {n_vars}/2 = {delta_Q:.4f}")
    print(f"  The radial part of the integral behaves like R^({d_L} + {d_C} + {d_Q} - 1 - p * {delta_Q:.2f}) = R^({power_Q} - p*{delta_Q:.2f})")
    print(f"  For convergence, {power_Q} - p*{delta_Q:.2f} < -1")
    print(f"  This implies p > ({power_Q} + 1) / {delta_Q:.4f} = {p_c_Q}")
    print("-" * 35)

    # Case 3: Linear terms dominate (Region D_L)
    # The integral behaves like integral(R^(d_Q+d_C) * (R^(-delta_L))^p * R^(d_L-1) dR)
    power_L = d_Q + d_C + d_L - 1
    # For convergence, power_L - p*delta_L < -1
    # p > (power_L + 1) / delta_L
    p_c_L = (power_L + 1) / delta_L
    print(f"In the region where linear terms dominate (D_L):")
    print(f"  Dimension of quadratic space d_Q = {d_Q}")
    print(f"  Dimension of cubic space d_C = {d_C}")
    print(f"  Dimension of linear space d_L = {d_L}")
    print(f"  Decay rate delta_L = {delta_L:.4f}")
    print(f"  The radial part of the integral behaves like R^({d_Q} + {d_C} + {d_L} - 1 - p * {delta_L:.2f}) = R^({power_L} - p*{delta_L:.2f})")
    print(f"  For convergence, {power_L} - p*{delta_L:.2f} < -1")
    print(f"  This implies p > ({power_L} + 1) / {delta_L:.4f} = {p_c_L}")
    print("-" * 35)

    # Step 4: The overall critical exponent is the maximum of the individual ones.
    p_c = max(p_c_C, p_c_Q, p_c_L)
    print(f"The overall integral converges if p is greater than the maximum of these thresholds.")
    print(f"p_c = max({p_c_C}, {p_c_Q}, {p_c_L}) = {p_c}")
    print("\nThe function I(a) is not in L^p for p <= p_c.")
    print(f"The largest p such that I is not in L^p(R^9) is {p_c}.")
    
    return p_c

if __name__ == '__main__':
    solve_critical_exponent()