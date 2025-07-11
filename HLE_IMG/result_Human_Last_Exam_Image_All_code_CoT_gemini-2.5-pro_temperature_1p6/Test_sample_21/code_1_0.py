def analyze_current_account_impact():
    """
    Analyzes the theoretical impact of a specific fiscal policy on a country's current account.

    The key identity is: ΔCA = ΔSp + ΔSg - ΔI
    where:
    ΔCA = Change in Current Account
    ΔSp = Change in Private Savings
    ΔSg = Change in Public Savings (ΔT - ΔG)
    ΔI  = Change in Investment
    """

    # 1. "a large fiscal expansion" means government spending (G) increases.
    # We use an illustrative value of 100 for the change in G (ΔG).
    delta_G = 100

    # 2. "which will be debt-financed" means taxes (T) do not change.
    delta_T = 0

    # 3. "increased government spending is completely offset by increased private saving".
    # This means the change in private savings (ΔSp) equals the change in government spending (ΔG).
    delta_Sp = delta_G

    # 4. We assume no change in investment (I) as it is not specified.
    delta_I = 0

    # First, calculate the change in public savings (ΔSg)
    delta_Sg = delta_T - delta_G

    # Now, calculate the total change in the current account (ΔCA)
    # ΔCA = ΔSp + ΔSg - ΔI
    delta_CA = delta_Sp + delta_Sg - delta_I

    # Print the explanation and the result
    print("The equation for the change in the Current Account (ΔCA) is:")
    print("ΔCA = ΔSp + (ΔT - ΔG) - ΔI")
    print("\nBased on the scenario:")
    print(f"- Change in Government Spending (ΔG): {delta_G}")
    print(f"- Change in Private Savings (ΔSp): {delta_Sp}")
    print(f"- Change in Taxes (ΔT): {delta_T}")
    print(f"- Change in Investment (ΔI): {delta_I}")

    print("\nPlugging the values into the equation:")
    # We explicitly show each number in the final calculation
    print(f"ΔCA = {delta_Sp} + ({delta_T} - {delta_G}) - {delta_I}")
    print(f"ΔCA = {delta_Sp} + ({delta_Sg}) - {delta_I}")
    print(f"ΔCA = {delta_CA}")

    print("\nConclusion: Under the given assumptions (full Ricardian equivalence),")
    print("the fiscal expansion has no net impact on the current account balance.")

analyze_current_account_impact()