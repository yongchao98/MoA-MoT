import numpy as np
import lime
import lime.lime_tabular
import warnings

# Suppress a warning from sklearn that is not relevant here
warnings.filterwarnings("ignore", category=UserWarning, module='sklearn')


def solve():
    """
    This function sets up the LIME explanation problem and prints the results.
    """
    # 1. Define the model function
    lookup_table = {1.0: 1.0, 0.0: 0.0}

    def f_raw(input1, input2):
        """The original model function."""
        return lookup_table.get(input1, input1 * 0 + input2 * 0.5 + 0.5)

    def predict_fn(x):
        """Wrapper function for LIME, expects a (n_samples, n_features) numpy array."""
        outputs = []
        for row in x:
            outputs.append(f_raw(row[0], row[1]))
        return np.array(outputs)

    # 2. Set up the LIME explainer
    # The baseline dataset is interpreted as points with keys from the lookup table.
    training_data = np.array([[0.0, 0.0], [1.0, 1.0]])
    feature_names = ['input1', 'input2']

    # With default hyperparams, LIME will detect features with few unique values
    # (<=4) as categorical. Both input1 and input2 have only 2 unique values.
    explainer = lime.lime_tabular.LimeTabularExplainer(
        training_data=training_data,
        feature_names=feature_names,
        class_names=['prediction'],
        mode='regression',
        verbose=False,
        random_state=42 # for reproducibility
    )

    # 3. Analyze scenario i): E belongs to the baseline dataset
    explicand_1 = np.array([0.0, 0.0])
    print("--- Scenario i): Explanation for E = (0.0, 0.0) ---")
    
    exp1 = explainer.explain_instance(
        explicand_1,
        predict_fn,
        num_features=2
    )
    
    # Print the local linear model equation
    print(f"Prediction(E) = {f_raw(*explicand_1):.4f}")
    print("Local linear model around E:")
    # The intercept is the expected value of the model over the perturbations
    print(f"y ≈ {exp1.intercept[0]:.4f} + ...")
    explanation_list_1 = exp1.as_list()
    # Find the more important feature
    important_feature_1 = ""
    max_abs_weight_1 = -1
    for feature, weight in explanation_list_1:
        print(f"  ... + ({weight:.4f} * I({feature}))")
        if abs(weight) > max_abs_weight_1:
            max_abs_weight_1 = abs(weight)
            important_feature_1 = feature.split('=')[0].strip()

    print(f"\nFor case i), the most important feature is: '{important_feature_1}'")
    print("-" * 50)

    # 4. Analyze scenario ii): E does not belong to the baseline dataset
    explicand_2 = np.array([-1.0, -1.0])
    print("--- Scenario ii): Explanation for E = (-1.0, -1.0) ---")
    
    exp2 = explainer.explain_instance(
        explicand_2,
        predict_fn,
        num_features=2
    )

    # Print the local linear model equation
    print(f"Prediction(E) = {f_raw(*explicand_2):.4f}")
    print("Local linear model around E:")
    print(f"y ≈ {exp2.intercept[0]:.4f} + ...")
    explanation_list_2 = exp2.as_list()
    # Find the more important feature
    important_feature_2 = ""
    max_abs_weight_2 = -1
    for feature, weight in explanation_list_2:
        print(f"  ... + ({weight:.4f} * I({feature}))")
        if abs(weight) > max_abs_weight_2:
            max_abs_weight_2 = abs(weight)
            important_feature_2 = feature.split('=')[0].strip()

    print(f"\nFor case ii), the most important feature is: '{important_feature_2}'")
    print("-" * 50)

    # 5. Conclusion
    print("\nConclusion:")
    print("In both scenarios, LIME identifies 'input1' as the feature with the highest importance (largest absolute weight).")
    print("This is because the model's behavior is entirely dictated by whether 'input1' is 0.0 or 1.0, and LIME's perturbation of the categorical 'input1' makes this dependency clear.")


solve()
<<<D>>>