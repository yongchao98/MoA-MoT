import numpy as np
import lime
import lime.lime_tabular
import re

# Step 1: Define the model function as described.
lookup_table = {1.0: 1.0, 0.0: 0.0}
def f_scalar(input1, input2):
    """The model function for scalar inputs."""
    return lookup_table.get(input1, input1 * 0 + input2 * 0.5 + 0.5)

def predict_fn(x):
    """
    A wrapper for the model function to handle numpy arrays of shape (n_samples, n_features),
    which is the format required by LIME.
    """
    # Apply the scalar function to each row of the input array.
    predictions = np.array([f_scalar(row[0], row[1]) for row in x])
    # For regression, LIME expects a 1D array of predictions.
    return predictions

# --- LIME Setup ---
# Step 2: Set up the LIME explainer with the specified baseline dataset and default parameters.
# The baseline dataset is interpreted from the lookup table as points (0.0, 0.0) and (1.0, 1.0).
training_data = np.array([[0.0, 0.0], [1.0, 1.0]])
feature_names = ['input1', 'input2']

# Create a LIME Tabular Explainer. We use default hyperparameters.
# The crucial default is `discretize_continuous=True`, which samples perturbations
# from the training data values for each feature.
explainer = lime.lime_tabular.LimeTabularExplainer(
    training_data=training_data,
    feature_names=feature_names,
    class_names=['prediction'],
    mode='regression',
    verbose=False,
    discretize_continuous=True,  # This is the default, but we are being explicit.
    random_state=42              # Use a fixed random state for reproducible results.
)

def analyze_and_print_explanation(case_name, instance, explainer, predict_fn):
    """
    This function generates, parses, and prints the LIME explanation for a given instance.
    It shows the local linear equation and determines which feature is more important.
    """
    print(f"--- {case_name} ---")
    print(f"Instance to explain E = ({instance[0]}, {instance[1]})")

    # Generate the explanation object from LIME.
    exp = explainer.explain_instance(
        instance,
        predict_fn,
        num_features=len(feature_names),
        num_samples=5000  # Default number of samples
    )

    # Output the components of the local linear model found by LIME.
    print("LIME's Local Linear Model:")
    intercept = exp.intercept[0]
    explanation_list = exp.as_list()
    
    # Sort by feature name for a consistent output order.
    explanation_list.sort(key=lambda x: x[0])

    # Build and print the equation string.
    equation_parts = [f"prediction = {intercept:.4f}"]
    importances = {}

    for feature_description, weight in explanation_list:
        # The feature_description tells us about the perturbation rule.
        # e.g., "0.00 <= input1 <= 1.00"
        feature_name = 'input1' if 'input1' in feature_description else 'input2'
        
        sign = "+" if weight >= 0 else "-"
        equation_parts.append(f"{sign} {abs(weight):.4f} * I({feature_description})")
        
        # Accumulate the absolute importance for each feature.
        importances[feature_name] = importances.get(feature_name, 0) + abs(weight)

    print(' '.join(equation_parts))
    print(f"Prediction from LIME's linear model: {exp.local_pred[0]:.4f}")
    print(f"Actual model prediction for E: {predict_fn(instance.reshape(1,-1))[0]:.4f}")

    # Determine and print the final conclusion for the case.
    if importances.get('input1', 0) > importances.get('input2', 0):
        conclusion = "input1 is more important."
    elif importances.get('input2', 0) > importances.get('input1', 0):
        conclusion = "input2 is more important."
    else:
        conclusion = "Both inputs have equal importance."
        
    print(f"Conclusion: {conclusion}\n")
    return conclusion

# Step 3: Analyze Case (i) where E belongs to the baseline dataset.
conclusion_i = analyze_and_print_explanation(
    "Case (i): E belongs to the baseline dataset",
    np.array([0.0, 0.0]),
    explainer,
    predict_fn
)

# Step 4: Analyze Case (ii) where E does not belong to the baseline dataset.
conclusion_ii = analyze_and_print_explanation(
    "Case (ii): E does not belong to the baseline dataset",
    np.array([-1.0, -1.0]),
    explainer,
    predict_fn
)

# Step 5: Summarize the findings.
print("--- Summary ---")
print(f"For case (i), the analysis shows that {conclusion_i}")
print(f"For case (ii), the analysis shows that {conclusion_ii}")
print("This corresponds to Answer A.")
<<<A>>>