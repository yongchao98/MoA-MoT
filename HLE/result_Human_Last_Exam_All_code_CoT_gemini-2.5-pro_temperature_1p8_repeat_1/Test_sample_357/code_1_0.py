def find_separating_equilibrium():
    """
    Calculates the pair of contracts in the separating equilibrium.
    """

    # Step 1: Determine the optimal value-added levels (v_L, v_H).
    # Based on the analysis, the principal chooses v_L to maximize the low-type's
    # social surplus (v_L - v_L^2), which occurs at v_L = 0.5.
    v_L = 0.5

    # The principal chooses v_H to minimize the information rent for the low type,
    # subject to constraints. This leads to v_H = 1.0.
    v_H = 1.0

    # Step 2: Determine the prices (p_L, p_H) based on binding constraints.
    # The high type's participation constraint is binding: p_H - v_H = 0
    p_H = v_H

    # The low type's incentive compatibility constraint is binding: p_L - v_L^2 = p_H - v_H^2
    p_L = v_L**2 + p_H - v_H**2

    # Step 3: Print the contracts, showing the calculation.
    print("The separating equilibrium consists of two contracts:")
    print("-" * 50)
    
    # Contract for the Low Type Handyman
    print("Contract for Low Type (L):")
    print(f"  Value added (v_L): {v_L}")
    # The prompt requests to output each number in the final equation.
    print(f"  Price (p_L) Calculation: p_L = v_L^2 + v_H - v_H^2 = {v_L}^2 + {v_H} - {v_H}^2")
    print(f"  Price (p_L): {p_L}")
    print(f"  So, the contract is (v_L, p_L) = ({v_L}, {p_L})")
    print("-" * 50)
    
    # Contract for the High Type Handyman
    print("Contract for High Type (H):")
    print(f"  Value added (v_H): {v_H}")
    print(f"  Price (p_H) Calculation: p_H = v_H = {v_H}")
    print(f"  Price (p_H): {p_H}")
    print(f"  So, the contract is (v_H, p_H) = ({v_H}, {p_H})")
    print("-" * 50)

if __name__ == '__main__':
    find_separating_equilibrium()
