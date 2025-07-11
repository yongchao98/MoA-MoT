import sympy as sp

def solve_contract_problem():
    """
    Calculates the separating equilibrium contracts for Chelsea and the handymen.
    """

    # Step 1 & 2: Define variables and constraints based on the analysis.
    # The problem is to maximize Chelsea's expected profit:
    # E[profit] = (5/6)*(v_H - p_H) + (1/6)*(v_L - p_L)
    # subject to IR and IC constraints.

    # Step 3: Solve the simplified optimization problem.
    # Our analysis showed that the IR constraints are binding, so:
    # p_L = v_L^2
    # p_H = v_H
    #
    # Substituting this into the expected profit function yields:
    # E[profit] = (1/6) * (v_L - v_L^2)
    #
    # The IC constraints under binding IRs become:
    # 0 <= v_L <= 1
    # v_H >= 1
    #
    # We need to maximize (v_L - v_L^2) with respect to v_L.
    v_L_var = sp.symbols('v_L')
    profit_L = v_L_var - v_L_var**2
    
    # Take the derivative with respect to v_L and solve for 0
    # to find the maximum.
    d_profit_L = sp.diff(profit_L, v_L_var)
    # 1 - 2*v_L = 0
    v_L_sol = sp.solve(d_profit_L, v_L_var)[0] # Result is 1/2

    # Step 4: Determine the contract parameters.
    v_L = float(v_L_sol)
    p_L = v_L**2
    
    # For the high type, v_H does not affect Chelsea's profit.
    # We choose the highest possible value according to the productivity constraint.
    # The valid range for v_H from IC and productivity constraints is [1, 30].
    v_H_max = 30
    v_H = float(v_H_max)
    p_H = v_H

    # Output the results
    print("The separating equilibrium is defined by the following pair of contracts:")
    print("\nContract for the Low-Type Handyman:")
    print(f"(v_L, p_L) = ({v_L}, {p_L})")
    print(f"Value added (v_L): {v_L}")
    print(f"Price (p_L): {p_L}")

    print("\nContract for the High-Type Handyman:")
    print(f"(v_H, p_H) = ({v_H}, {p_H})")
    print(f"Value added (v_H): {v_H}")
    print(f"Price (p_H): {p_H}")

solve_contract_problem()