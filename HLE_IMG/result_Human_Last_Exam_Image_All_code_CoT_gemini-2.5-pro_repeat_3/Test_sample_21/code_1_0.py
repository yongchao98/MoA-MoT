def analyze_current_account():
    """
    Analyzes the impact of a fiscal expansion on the current account
    based on the given assumptions.
    """
    print("This analysis uses the national saving and investment identity:")
    print("Current Account (CA) = National Saving (S) - Investment (I)")
    print("Where National Saving (S) = Private Saving (Sp) + Public Saving (Sg)\n")

    # Let's represent the large fiscal expansion with a symbolic value, e.g., 100 units.
    fiscal_expansion_delta_G = 100

    print(f"Step 1: A debt-financed fiscal expansion of {fiscal_expansion_delta_G} units occurs.")
    print("This means Government Spending (G) increases by 100, and Taxes (T) are unchanged.")
    print("Public Saving (Sg = T - G) decreases by the same amount.")
    delta_Sg = -fiscal_expansion_delta_G
    print(f"Change in Public Saving (ΔSg) = {delta_Sg}\n")

    print(f"Step 2: The problem assumes this is completely offset by increased private saving.")
    delta_Sp = fiscal_expansion_delta_G
    print(f"Change in Private Saving (ΔSp) = {delta_Sp}\n")

    print("Step 3: Calculate the total change in National Saving (ΔS).")
    print("ΔS = ΔSp + ΔSg")
    delta_S_national = delta_Sp + delta_Sg
    print(f"Change in National Saving (ΔS) = {delta_Sp} + ({delta_Sg}) = {delta_S_national}\n")

    print("Step 4: Assume there is no change in Investment (I).")
    delta_I = 0
    print(f"Change in Investment (ΔI) = {delta_I}\n")

    print("Step 5: Calculate the final change in the Current Account (ΔCA).")
    print("ΔCA = ΔS - ΔI")
    delta_CA = delta_S_national - delta_I
    
    # Final output showing the equation with numbers
    print("\n--- Final Calculation ---")
    print("The change in the current account is given by the equation:")
    print(f"ΔCA = (Change in Private Saving + Change in Public Saving) - Change in Investment")
    print(f"ΔCA = ({delta_Sp} + ({delta_Sg})) - {delta_I}")
    print(f"ΔCA = {delta_S_national} - {delta_I}")
    print(f"ΔCA = {delta_CA}")

    print("\nConclusion: There is no impact on the country's current account balance.")

analyze_current_account()