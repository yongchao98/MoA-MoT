import math

def solve_contracts():
    """
    This function calculates and prints the separating equilibrium contracts.
    The derivation is explained in the text above.
    """
    # Based on the derivation, a possible equilibrium is found by choosing v_L and v_H
    # that satisfy a set of constraints. We select a simple integer solution.
    v_L = 1.0
    v_H = 0.0

    # Prices are determined by the binding IR-H and IC-H constraints
    # p_H = v_H + 30
    # p_L = v_L + 30
    p_L = v_L + 30.0
    p_H = v_H + 30.0

    print("A separating equilibrium is constituted by the following pair of contracts:")
    print("Contract for Low Type (v_L, p_L):")
    print(f"Value added (v_L): {v_L}")
    print(f"Price (p_L): {p_L}")
    print("\nContract for High Type (v_H, p_H):")
    print(f"Value added (v_H): {v_H}")
    print(f"Price (p_H): {p_H}")

    # Final answer formatting
    # The final output requires printing the equations for the contracts.
    # So we format it as requested.
    v_L_int = int(v_L)
    p_L_int = int(p_L)
    v_H_int = int(v_H)
    p_H_int = int(p_H)
    
    answer = f"(v_L, p_L) = ({v_L_int}, {p_L_int}), (v_H, p_H) = ({v_H_int}, {p_H_int})"
    print(f"\nFinal Answer in the required format: <<<({v_L_int}, {p_L_int}), ({v_H_int}, {p_H_int})>>>")


solve_contracts()