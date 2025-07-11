import numpy as np
import lime
import lime.lime_tabular
import re

def solve():
    """
    This script demonstrates LIME's behavior on a given function for two scenarios,
    determining which feature is more important in each case.
    """
    # Define the model function as described in the problem.
    lookup_table = {1.0: 1.0, 0.0: 0.0}
    def f(input1, input2):
        # Use float() to ensure numpy floats match the dictionary keys
        input1_float = float(input1)
        if input1_float in lookup_table:
            return lookup_table[input1_float]
        else:
            return input1 * 0 + input2 * 0.5 + 0.5

    # Create a prediction function wrapper that LIME can use.
    # It must accept a numpy array (n_samples, n_features) and return an array (n_samples,).
    def predict_fn(x):
        return np.array([f(row[0], row[1]) for row in x])

    # The baseline dataset consists of points where features are 0.0 or 1.0.
    training_data = np.array([[0.0, 0.0], [0.0, 1.0], [1.0, 0.0], [1.0, 1.0]], dtype=np.float64)
    feature_names = ['input1', 'input2']

    # Initialize the LIME explainer. We use continuous perturbations to analyze
    # the function's local behavior, which is central to LIME's purpose.
    explainer = lime.lime_tabular.LimeTabularExplainer(
        training_data=training_data,
        feature_names=feature_names,
        mode='regression',
        discretize_continuous=False,  # Use continuous perturbations
        random_state=42  # For reproducible results
    )

    # --- Scenario i) E belongs to the baseline dataset: E = (0.0, 0.0) ---
    print("--- Case i: Explaining for E = (0.0, 0.0) ---")
    explicand1 = np.array([0.0, 0.0])
    exp1 = explainer.explain_instance(
        explicand1,
        predict_fn,
        num_features=2,
        num_samples=5000
    )
    
    print("LIME Explanation (Local Linear Model):")
    # Output each term of the explanation equation
    explanation_list1 = exp1.as_list()
    for feat_explanation, weight in explanation_list1:
        print(f"Term: ({feat_explanation}) * Weight: ({weight:.4f})")

    # Determine the more important feature by comparing absolute weights
    weights1 = dict(exp1.local_exp[1])
    importance1_input1 = abs(weights1.get(0, 0.0))
    importance1_input2 = abs(weights1.get(1, 0.0))
    conclusion1 = "input1" if importance1_input1 > importance1_input2 else "input2"
    print(f"\nConclusion for Case i: '{conclusion1}' is more important.")
    print("="*50 + "\n")


    # --- Scenario ii) E does not belong to the baseline dataset: E = (-1.0, -1.0) ---
    print("--- Case ii: Explaining for E = (-1.0, -1.0) ---")
    explicand2 = np.array([-1.0, -1.0])
    exp2 = explainer.explain_instance(
        explicand2,
        predict_fn,
        num_features=2,
        num_samples=5000
    )

    print("LIME Explanation (Local Linear Model):")
    # Output each term of the explanation equation
    explanation_list2 = exp2.as_list()
    for feat_explanation, weight in explanation_list2:
        print(f"Term: ({feat_explanation}) * Weight: ({weight:.4f})")
    
    # Determine the more important feature by comparing absolute weights
    weights2 = dict(exp2.local_exp[1])
    importance2_input1 = abs(weights2.get(0, 0.0))
    importance2_input2 = abs(weights2.get(1, 0.0))
    conclusion2 = "input1" if importance2_input1 > importance2_input2 else "input2"
    print(f"\nConclusion for Case ii: '{conclusion2}' is more important.")

solve()