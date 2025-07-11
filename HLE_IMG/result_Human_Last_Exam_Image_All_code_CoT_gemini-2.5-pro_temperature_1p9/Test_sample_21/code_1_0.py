def analyze_current_account():
    """
    Analyzes the impact of a debt-financed fiscal expansion on the current account
    under the condition that private saving completely offsets government spending.
    """

    # 1. Define the initial change from the fiscal expansion.
    # Let's use a hypothetical value of 100 units for the increase in government spending (G).
    delta_G = 100

    # 2. Determine the impact on public saving (Sg).
    # Public saving = Taxes (T) - Government Spending (G).
    # Since the expansion is debt-financed, T is unchanged.
    # The change in public saving is the negative of the change in G.
    delta_Sg = -delta_G

    # 3. Determine the impact on private saving (Sp).
    # The problem states increased private saving completely offsets the increased government spending.
    delta_Sp = delta_G

    # 4. Calculate the total change in National Saving (S).
    # National Saving (S) = Private Saving (Sp) + Public Saving (Sg).
    # The change in S (delta_S) is the sum of the changes in Sp and Sg.
    delta_S = delta_Sp + delta_Sg

    # 5. Assume no change in investment (I), as is typical in this theoretical scenario.
    delta_I = 0

    # 6. Calculate the final change in the Current Account Balance (NX).
    # The identity is: Current Account (NX) = National Saving (S) - Investment (I).
    # The change in the current account (delta_NX) is calculated from the changes in S and I.
    delta_NX = delta_S - delta_I

    # --- Print the step-by-step analysis ---
    print("--- Theoretical Impact of Fiscal Expansion on Current Account ---")
    print("\nThe core relationship is: ΔNX = ΔS - ΔI")
    print("where ΔNX is change in Current Account, ΔS is change in National Saving, and ΔI is change in Investment.")
    print("\nNational Saving (S) = Private Saving (Sp) + Public Saving (Sg).")
    print("-" * 60)
    
    print(f"Step 1: Fiscal expansion implies an increase in government spending (ΔG) of {delta_G}.")
    print(f"Step 2: Debt-financing this reduces public saving (ΔSg) by {abs(delta_Sg)}, so ΔSg = {delta_Sg}.")
    print(f"Step 3: Private saving increases to offset this (ΔSp), so ΔSp = {delta_Sp}.")
    print("-" * 60)
    
    print("Step 4: Calculate the total change in National Saving (ΔS):")
    print(f"   ΔS = ΔSp + ΔSg")
    print(f"   ΔS = {delta_Sp} + ({delta_Sg})")
    print(f"   ΔS = {delta_S}")
    print("-" * 60)
    
    print("Step 5: Calculate the final change in the Current Account (ΔNX):")
    print(f"   We assume no change in investment, so ΔI = {delta_I}.")
    print(f"   ΔNX = ΔS - ΔI")
    print(f"   ΔNX = {delta_S} - {delta_I}")
    print(f"   ΔNX = {delta_NX}")
    print("-" * 60)

    print("\nConclusion: The total change in the current account is zero.")


analyze_current_account()