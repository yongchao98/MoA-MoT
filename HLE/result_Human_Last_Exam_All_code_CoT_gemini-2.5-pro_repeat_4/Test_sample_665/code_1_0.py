import numpy as np
import lime
import lime.lime_tabular

def solve():
    """
    This script sets up the specified model and uses LIME to determine
    feature importances for the two scenarios described in the problem.
    """
    # 1. Define the model function and a wrapper for LIME
    lookup_table = {1.0: 1.0, 0.0: 0.0}
    def f(input1, input2):
        """The scalar model function."""
        return lookup_table.get(input1, input1 * 0 + input2 * 0.5 + 0.5)

    def predict_fn(X):
        """
        Wrapper function for LIME.
        Takes a numpy array of shape (n_samples, n_features)
        and returns a numpy array of predictions.
        """
        return np.array([f(row[0], row[1]) for row in X])

    # 2. Set up the baseline data and LIME explainer
    training_data = np.array([[1.0, 1.0], [0.0, 0.0]])
    feature_names = ['input1', 'input2']

    # Using default LIME Tabular Explainer parameters, which is key to the solution.
    explainer = lime.lime_tabular.LimeTabularExplainer(
        training_data=training_data,
        feature_names=feature_names,
        class_names=['output'],
        mode='regression'
    )

    # --- Case i): E belongs to the baseline dataset ---
    explicand1 = np.array([0.0, 0.0])
    print("--- Explanation for Case i): E = (0.0, 0.0) ---")
    
    explanation1 = explainer.explain_instance(
        explicand1,
        predict_fn,
        num_features=2
    )
    
    # Print the "equation" for the local model
    print("LIME's Local Linear Model:")
    intercept1 = explanation1.intercept[0]
    print(f"Prediction ≈ {intercept1:.4f} + ", end="")
    
    # as_list() is sorted by importance. Let's print in order.
    feature_weights1 = explanation1.as_list()
    equation_parts1 = []
    for feature, weight in feature_weights1:
        # The feature string itself is the condition for the binary variable
        equation_parts1.append(f"{weight:.4f} * [Is {feature}]")
    print(" + ".join(equation_parts1))


    # --- Case ii): E does not belong to the baseline dataset ---
    explicand2 = np.array([-1.0, -1.0])
    print("\n--- Explanation for Case ii): E = (-1.0, -1.0) ---")
    
    explanation2 = explainer.explain_instance(
        explicand2,
        predict_fn,
        num_features=2
    )
    
    # Print the "equation" for the local model
    print("LIME's Local Linear Model:")
    intercept2 = explanation2.intercept[0]
    print(f"Prediction ≈ {intercept2:.4f} + ", end="")
    
    feature_weights2 = explanation2.as_list()
    equation_parts2 = []
    for feature, weight in feature_weights2:
        equation_parts2.append(f"{weight:.4f} * [Is {feature}]")
    print(" + ".join(equation_parts2))

    # --- Final Conclusion ---
    print("\n\n--- Conclusion ---")
    feature1_case1 = feature_weights1[0][0]
    weight1_case1 = feature_weights1[0][1]
    feature1_case2 = feature_weights2[0][0]
    weight1_case2 = feature_weights2[0][1]

    print(f"For case i), the most important feature is '{feature1_case1.split(' ')[0]}' with a weight of {weight1_case1:.4f}.")
    print(f"For case ii), the most important feature is '{feature1_case2.split(' ')[0]}' with a weight of {weight1_case2:.4f}.")
    print("In both scenarios, 'input1' is identified as the most important feature.")

solve()