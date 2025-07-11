def analyze_fiscal_impact():
    """
    Analyzes the impact of a debt-financed fiscal expansion on the current account
    under the assumption of a full private saving offset (Ricardian Equivalence).
    """
    # Let's use a hypothetical value for the increase in government spending to illustrate.
    delta_G = 100

    print("Step 1: A large fiscal expansion occurs.")
    print(f"Change in Government Spending (ΔG) = +{delta_G}\n")

    # The expansion is debt-financed, so taxes (T) do not change.
    delta_T = 0
    # The change in Public Saving (Sg = T - G) is ΔSg = ΔT - ΔG.
    delta_Sg = delta_T - delta_G
    print("Step 2: The expansion is debt-financed, so Public Saving (Sg) decreases.")
    print(f"Change in Public Saving (ΔSg) = {delta_T} - {delta_G} = {delta_Sg}\n")

    # The core assumption: increased private saving completely offsets the increased government spending.
    # This means the change in Private Saving (ΔSp) is equal to the change in Government Spending (ΔG).
    delta_Sp = delta_G
    print("Step 3: Private Saving (Sp) increases to offset the fiscal expansion.")
    print(f"Change in Private Saving (ΔSp) = +{delta_Sp}\n")

    # National Saving (S) is the sum of private and public saving (S = Sp + Sg).
    # The total change in National Saving (ΔS) is the sum of the changes in its components.
    delta_S = delta_Sp + delta_Sg
    print("Step 4: Calculate the total change in National Saving (S = Sp + Sg).")
    print(f"Change in National Saving (ΔS) = ΔSp + ΔSg = {delta_Sp} + ({delta_Sg}) = {delta_S}\n")

    # The Current Account balance (CA) is National Saving (S) minus Investment (I).
    # We assume Investment (I) is unchanged, as no information suggests otherwise.
    delta_I = 0
    delta_CA = delta_S - delta_I
    print("Step 5: Calculate the final impact on the Current Account (CA = S - I).")
    print("Assuming Investment (I) is unchanged (ΔI = 0).")
    print(f"Change in Current Account (ΔCA) = ΔS - ΔI = {delta_S} - {delta_I} = {delta_CA}\n")

    print("Conclusion: The decrease in public saving is perfectly matched by an increase in private saving.")
    print("This results in no change to national saving, and therefore no change to the current account balance.")

analyze_fiscal_impact()