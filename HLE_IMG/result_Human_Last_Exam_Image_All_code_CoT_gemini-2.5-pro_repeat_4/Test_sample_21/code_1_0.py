def analyze_current_account():
    """
    Analyzes the impact of a specific fiscal expansion on the current account balance.
    """
    # 1. Assume a value for the fiscal expansion (increase in Government Spending, G).
    # Let's say the increase in G is 100 units.
    delta_G = 100
    print(f"Step 1: A debt-financed fiscal expansion increases government spending (G).")
    print(f"Let's assume the change in G (ΔG) is: {delta_G}\n")

    # 2. Calculate the change in government saving (Sg).
    # Government saving (Sg) = Taxes (T) - Government Spending (G).
    # Since G increases and T is constant, Sg decreases by the same amount.
    delta_Sg = -delta_G
    print(f"Step 2: This reduces government saving (Sg). The change in Sg (ΔSg) is the negative of ΔG.")
    print(f"ΔSg = {delta_Sg}\n")

    # 3. Determine the change in private saving (Sp).
    # The problem states that increased private saving completely offsets the increased government spending.
    delta_Sp = delta_G
    print(f"Step 3: Private saving (Sp) increases to offset the fiscal expansion.")
    print(f"The change in Sp (ΔSp) is: {delta_Sp}\n")

    # 4. Calculate the change in national saving (S).
    # National Saving (S) = Private Saving (Sp) + Government Saving (Sg).
    # The change in national saving (ΔS) is the sum of the changes in Sp and Sg.
    delta_S = delta_Sp + delta_Sg
    print(f"Step 4: Calculate the total change in National Saving (ΔS).")
    print(f"ΔS = ΔSp + ΔSg")
    print(f"ΔS = {delta_Sp} + ({delta_Sg})")
    print(f"ΔS = {delta_S}\n")

    # 5. Determine the change in investment (I).
    # Under full Ricardian Equivalence, the policy does not affect interest rates, so investment is unchanged.
    delta_I = 0
    print(f"Step 5: Investment (I) is assumed to be unaffected by this policy.")
    print(f"The change in I (ΔI) is: {delta_I}\n")

    # 6. Calculate the final change in the Current Account (CA).
    # Current Account (CA) = National Saving (S) - Investment (I).
    # The change in the current account (ΔCA) is ΔS - ΔI.
    delta_CA = delta_S - delta_I
    print(f"Step 6: Calculate the final change in the Current Account (ΔCA).")
    print(f"ΔCA = ΔS - ΔI")
    print(f"ΔCA = {delta_S} - {delta_I}")
    print(f"ΔCA = {delta_CA}\n")

    print("Conclusion: The decrease in government saving is perfectly cancelled out by the increase in private saving.")
    print("As a result, national saving does not change, and the current account balance is unaffected.")


analyze_current_account()