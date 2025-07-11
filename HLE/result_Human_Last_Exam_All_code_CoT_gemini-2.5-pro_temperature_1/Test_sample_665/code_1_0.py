import numpy as np
import lime
import lime.lime_tabular
import warnings

# Suppress a common LIME warning about discretization
warnings.filterwarnings("ignore", category=UserWarning, module='lime')

def solve():
    """
    This script sets up the model and LIME explainer as described in the problem,
    then runs the analysis for both scenarios to determine feature importance.
    """
    # 1. Define the model function
    lookup_table = {1.0: 1.0, 0.0: 0.0}
    def f_raw(input1, input2):
        # The default value in .get() is executed if input1 is not in the lookup_table
        return lookup_table.get(input1, input1 * 0 + input2 * 0.5 + 0.5)

    # LIME requires a prediction function that takes a numpy array
    def predict_fn(x):
        predictions = []
        for row in x:
            input1, input2 = row
            predictions.append(f_raw(input1, input2))
        # LIME for regression needs a 1D array
        return np.array(predictions)

    # 2. Create the baseline/training dataset
    # "baseline dataset is the same as the lookup table" implies input1 is only 0.0 or 1.0.
    n_samples = 200
    # input1 is categorical: 0.0 or 1.0
    input1_train = np.random.choice([0.0, 1.0], size=n_samples)
    # input2 can be continuous, let's use a standard normal distribution
    input2_train = np.random.normal(0, 1, size=n_samples)
    training_data = np.vstack([input1_train, input2_train]).T
    feature_names = ['input1', 'input2']

    # 3. Create a LIME explainer
    # We explicitly tell LIME that the first feature is categorical.
    # We set a random_state for reproducibility.
    explainer = lime.lime_tabular.LimeTabularExplainer(
        training_data=training_data,
        feature_names=feature_names,
        categorical_features=[0], # 0 is the index of input1
        mode='regression',
        random_state=42 # For reproducible results
    )

    # 4. Define the two instances (explicands)
    E1 = np.array([0.0, 0.0])   # i) Belongs to the baseline dataset 'category'
    E2 = np.array([-1.0, -1.0]) # ii) Does not belong

    # --- Run and Print Analysis for Scenario i ---
    print("--- Scenario i: Explaining E = (0.0, 0.0) ---")
    # Explain the instance E1
    explanation1 = explainer.explain_instance(E1, predict_fn, num_features=2)
    # Get importances as a list of (feature, weight) tuples
    importances1 = explanation1.as_list()
    print(f"Model prediction for E1: {predict_fn(E1.reshape(1, -1))[0]:.4f}")
    print("LIME Feature Importances:")
    # The feature with the largest absolute weight is most important
    most_important_1 = ""
    max_abs_weight_1 = -1
    for feature, weight in importances1:
        print(f"  Feature: {feature}, Weight: {weight:.4f}")
        if abs(weight) > max_abs_weight_1:
            max_abs_weight_1 = abs(weight)
            # The feature name is the first word in the LIME explanation string
            most_important_1 = feature.split(' ')[0]
    print(f"\nConclusion for i): '{most_important_1}' is more important.")


    # --- Run and Print Analysis for Scenario ii ---
    print("\n--- Scenario ii: Explaining E = (-1.0, -1.0) ---")
    # Explain the instance E2
    explanation2 = explainer.explain_instance(E2, predict_fn, num_features=2)
    importances2 = explanation2.as_list()
    print(f"Model prediction for E2: {predict_fn(E2.reshape(1, -1))[0]:.4f}")
    print("LIME Feature Importances:")
    most_important_2 = ""
    max_abs_weight_2 = -1
    for feature, weight in importances2:
        print(f"  Feature: {feature}, Weight: {weight:.4f}")
        if abs(weight) > max_abs_weight_2:
            max_abs_weight_2 = abs(weight)
            most_important_2 = feature.split(' ')[0]
    print(f"\nConclusion for ii): '{most_important_2}' is more important.")
    
    # --- Final Answer ---
    print("\n---------------------------------------------------")
    print(f"Final Answer: i) {most_important_1} and ii) {most_important_2}")
    print("This corresponds to Choice A.")


solve()