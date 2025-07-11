import sys

# Suppress writing bytecode files
sys.dont_write_bytecode = True

def analyze_current_account():
    """
    Analyzes the impact of a specific fiscal policy on the current account balance.
    """
    print("To determine the impact, we use the national income identity for an open economy:")
    print("Current Account (CA) = National Savings (S) - Investment (I)")
    print("\nNational Savings can be split into Private Savings (S_p) and Public Savings (S_g):")
    print("CA = S_p + S_g - I")
    print("\nWe analyze the change (Δ) in each component based on the given scenario:")
    print("ΔCA = ΔS_p + ΔS_g - ΔI")
    print("------------------------------------------------------------------")

    # Let's use a hypothetical value for the fiscal expansion to demonstrate.
    # Assume the increase in government spending (ΔG) is 100 units.
    increase_in_G = 100

    # 1. Change in Public Savings (ΔS_g)
    # The fiscal expansion is debt-financed, meaning the government deficit increases.
    # Public savings (S_g = Taxes - G) decrease by the amount of the spending increase.
    delta_S_g = -increase_in_G
    print(f"1. A debt-financed fiscal expansion increases government spending by {increase_in_G}.")
    print(f"   This causes Public Savings (S_g) to decrease by the same amount.")
    print(f"   => ΔS_g = {delta_S_g}")

    # 2. Change in Private Savings (ΔS_p)
    # The problem states that increased private saving completely offsets the increased government spending.
    # This is the key assumption of Ricardian Equivalence.
    delta_S_p = increase_in_G
    print(f"\n2. Private Savings (S_p) increase to fully offset the increase in government spending.")
    print(f"   => ΔS_p = {delta_S_p}")

    # 3. Change in Investment (ΔI)
    # The problem provides no information about a change in investment.
    delta_I = 0
    print(f"\n3. We assume there is no change in Investment (I).")
    print(f"   => ΔI = {delta_I}")
    print("------------------------------------------------------------------")

    # Calculate the total change in the Current Account (ΔCA)
    delta_CA = delta_S_p + delta_S_g - delta_I

    print("\nNow, we calculate the total change in the Current Account using the equation:")
    print(f"ΔCA = ΔS_p + ΔS_g - ΔI")
    # Using f-string to explicitly show the values in the equation
    print(f"ΔCA = ({delta_S_p}) + ({delta_S_g}) - ({delta_I})")
    print(f"ΔCA = {delta_CA}")

    print("\nConclusion:")
    print("The increase in private savings exactly cancels out the decrease in public savings.")
    print("Therefore, national savings remain unchanged. With no change in investment, the fiscal")
    print("expansion has no impact on the country's current account balance.")


if __name__ == '__main__':
    analyze_current_account()