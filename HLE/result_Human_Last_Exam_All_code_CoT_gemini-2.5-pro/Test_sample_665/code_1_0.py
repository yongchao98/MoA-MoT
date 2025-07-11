import numpy as np
import lime
import lime.lime_tabular
import re

# Define the model function as specified
lookup_table = {1.0: 1.0, 0.0: 0.0}
def f(input1, input2):
    """The model function to be explained."""
    return lookup_table.get(input1, input1 * 0 + input2 * 0.5 + 0.5)

# Create a prediction function compatible with LIME, which expects a NumPy array
def predict_fn(x):
    """Wrapper for the model function to handle batch predictions."""
    outputs = []
    for row in x:
        outputs.append(f(row[0], row[1]))
    return np.array(outputs)

# Define the baseline dataset (used as training data for LIME)
# As per the problem, it's related to the lookup table keys
training_data = np.array([[0.0, 0.0], [1.0, 1.0], [0.0, 1.0], [1.0, 0.0]])
feature_names = ['input1', 'input2']

# Initialize the LIME Tabular Explainer with default hyperparameters
explainer = lime.lime_tabular.LimeTabularExplainer(
    training_data=training_data,
    feature_names=feature_names,
    mode='regression',
    random_state=42  # for reproducibility
)

def get_feature_importances(explanation):
    """Parses the explanation to get feature weights."""
    importances = {'input1': 0.0, 'input2': 0.0}
    for feature_condition, weight in explanation.as_list():
        if 'input1' in feature_condition:
            importances['input1'] += abs(weight)
        elif 'input2' in feature_condition:
            importances['input2'] += abs(weight)
    return importances

def print_explanation_summary(case_name, instance, explanation):
    """Prints a formatted summary of the LIME explanation."""
    print(f"--- {case_name} ---")
    print(f"Instance to explain E = {instance}")
    original_prediction = predict_fn(instance.reshape(1, -1))[0]
    print(f"Model prediction for E: f({instance[0]}, {instance[1]}) = {original_prediction:.4f}")
    
    print("\nLocal linear model (the 'equation'):")
    print(f"Intercept (base value) = {explanation.intercept[0]:.4f}")
    print("Feature weights:")
    for feature, weight in explanation.as_list():
        print(f"  '{feature}': {weight:.4f}")

    importances = get_feature_importances(explanation)
    
    print("\nConclusion:")
    if importances['input1'] > importances['input2']:
        print(f"input1 is more important (abs weight: {importances['input1']:.4f}) than input2 (abs weight: {importances['input2']:.4f})")
    elif importances['input2'] > importances['input1']:
        print(f"input2 is more important (abs weight: {importances['input2']:.4f}) than input1 (abs weight: {importances['input1']:.4f})")
    else:
        print("Both inputs have equal importance.")
    print("-" * (len(case_name) + 8))


# --- Case i): E belongs to the baseline dataset ---
instance_in_baseline = np.array([0.0, 0.0])
exp_in = explainer.explain_instance(
    instance_in_baseline,
    predict_fn,
    num_features=2
)
print_explanation_summary("Case i): E belongs to the baseline dataset", instance_in_baseline, exp_in)

print("\n\n")

# --- Case ii): E does not belong to the baseline dataset ---
instance_out_baseline = np.array([-1.0, -1.0])
exp_out = explainer.explain_instance(
    instance_out_baseline,
    predict_fn,
    num_features=2
)
print_explanation_summary("Case ii): E does not belong to the baseline dataset", instance_out_baseline, exp_out)