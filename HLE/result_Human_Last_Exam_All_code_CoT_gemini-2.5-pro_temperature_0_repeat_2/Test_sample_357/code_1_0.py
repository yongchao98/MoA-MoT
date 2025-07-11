import sympy

def solve_contracts():
    """
    This function solves for the separating equilibrium contracts in the given
    handyman problem.
    """
    # Define symbols for the variables
    v_L, v_H = sympy.symbols('v_L v_H')

    # --- Step 1: Maximize the profit function for the low type ---
    # The principal's profit from the low type, after substituting p_L = v_L**2, is v_L - v_L**2.
    # We find the maximum of this function.
    profit_L = v_L - v_L**2
    # Differentiate with respect to v_L and solve for 0 to find the maximum.
    d_profit_L = sympy.diff(profit_L, v_L)
    optimal_v_L = sympy.solve(d_profit_L, v_L)[0]

    # --- Step 2: Maximize the profit function for the high type ---
    # The principal's profit from the high type, after substituting p_H = v_H**2, is v_H - v_H**2.
    # This must be maximized subject to the constraint IR_H: v_H^2 - v_H >= 0, which implies v_H >= 1.
    # The function v_H - v_H**2 is decreasing for v_H > 0.5.
    # Therefore, to maximize it for v_H >= 1, we must choose the smallest possible value.
    optimal_v_H = 1

    # --- Step 3: Calculate the prices based on the optimal v values ---
    # p_L = v_L^2
    optimal_p_L = optimal_v_L**2
    # p_H = v_H^2
    optimal_p_H = optimal_v_H**2

    # --- Step 4: Print the results ---
    print("The separating equilibrium consists of two contracts: (v_L, p_L) and (v_H, p_H).")
    print("\nContract for the Low-Type Handyman:")
    print(f"Value added (v_L): {optimal_v_L}")
    print(f"Price (p_L): {optimal_p_L}")
    
    print("\nContract for the High-Type Handyman:")
    print(f"Value added (v_H): {optimal_v_H}")
    print(f"Price (p_H): {optimal_p_H}")

    # Final answer format
    # The pair of contracts is ((v_L, p_L), (v_H, p_H))
    # ((0.5, 0.25), (1, 1))
    
solve_contracts()