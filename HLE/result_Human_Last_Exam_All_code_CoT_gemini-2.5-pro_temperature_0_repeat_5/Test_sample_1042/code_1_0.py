import numpy as np

def demonstrate_non_identifiability():
    """
    This function demonstrates the non-identifiability of E(Y^a | A, L)
    by constructing two causal models that are indistinguishable from
    observational data but yield different counterfactual predictions.

    Let A, L, U be binary variables {0, 1}.
    Let Y be a continuous outcome.
    Let epsilon be a standard normal noise term N(0,1).

    In both models, we assume the unmeasured confounder U causes both A and Y.
    """

    # --- Define Model 1 ---
    # In this model, the unmeasured confounder U has a positive effect on Y.
    # Y^a = 2*L + a + U
    def model_1_counterfactual_expectation(a, A_obs, L_obs):
        # The causal relationship is A = U.
        # Therefore, E[U | A_obs] = A_obs.
        # E[Y^a | A_obs, L_obs] = E[2*L_obs + a + U + epsilon | A_obs, L_obs]
        #                       = 2*L_obs + a + E[U | A_obs]
        #                       = 2*L_obs + a + A_obs
        return 2 * L_obs + a + A_obs

    # --- Define Model 2 ---
    # In this model, the unmeasured confounder U has a negative effect on Y.
    # Y^a = 2*L + 3*a - U
    def model_2_counterfactual_expectation(a, A_obs, L_obs):
        # The causal relationship is A = U.
        # Therefore, E[U | A_obs] = A_obs.
        # E[Y^a | A_obs, L_obs] = E[2*L_obs + 3*a - U + epsilon | A_obs, L_obs]
        #                       = 2*L_obs + 3*a - E[U | A_obs]
        #                       = 2*L_obs + 3*a - A_obs
        return 2 * L_obs + 3 * a - A_obs

    # --- Check Observational Equivalence ---
    # The observed outcome Y occurs when the intervention 'a' is the same as
    # the observed treatment 'A'. So we set a = A.
    # For Model 1: E[Y | A, L] = 2*L + A + A = 2*L + 2*A
    # For Model 2: E[Y | A, L] = 2*L + 3*A - A = 2*L + 2*A
    # Both models produce the same expectation for the observed data.

    # Let's pick some specific values for observed variables
    A_obs = 0
    L_obs = 1

    # Calculate the observed expectation E(Y | A, L) for both models
    obs_exp_m1 = model_1_counterfactual_expectation(a=A_obs, A_obs=A_obs, L_obs=L_obs)
    obs_exp_m2 = model_2_counterfactual_expectation(a=A_obs, A_obs=A_obs, L_obs=L_obs)

    print("--- Demonstrating Non-Identifiability ---")
    print(f"Let's analyze the subpopulation where observed treatment A = {A_obs} and observed confounder L = {L_obs}.\n")

    print("1. Check for Observational Equivalence:")
    print(f"   The expected observed outcome E(Y | A={A_obs}, L={L_obs}) in Model 1 is: {obs_exp_m1}")
    print(f"   The expected observed outcome E(Y | A={A_obs}, L={L_obs}) in Model 2 is: {obs_exp_m2}")
    print("   The models are observationally indistinguishable.\n")

    # --- Check Counterfactual Difference ---
    # Now, let's ask a counterfactual question: What would the outcome be for this
    # same group if their treatment had been a=1 instead of A=0?
    a_cf = 1

    cf_exp_m1 = model_1_counterfactual_expectation(a=a_cf, A_obs=A_obs, L_obs=L_obs)
    cf_exp_m2 = model_2_counterfactual_expectation(a=a_cf, A_obs=A_obs, L_obs=L_obs)

    print("2. Check for Counterfactual Difference:")
    print(f"   We want to identify E(Y^a | A, L) for a = {a_cf}, A = {A_obs}, L = {L_obs}.")
    print(f"   Model 1 predicts: E(Y^a={a_cf} | A={A_obs}, L={L_obs}) = 2 * {L_obs} + {a_cf} + {A_obs} = {cf_exp_m1}")
    print(f"   Model 2 predicts: E(Y^a={a_cf} | A={A_obs}, L={L_obs}) = 2 * {L_obs} + 3 * {a_cf} - {A_obs} = {cf_exp_m2}")
    print("\nSince the same observed data can lead to different answers for the counterfactual quantity,")
    print("the quantity E(Y^a | A, L) is not identifiable.")

if __name__ == '__main__':
    demonstrate_non_identifiability()