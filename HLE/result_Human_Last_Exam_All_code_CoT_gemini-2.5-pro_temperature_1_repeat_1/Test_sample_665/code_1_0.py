import numpy as np
import lime
import lime.lime_tabular
import warnings

# Suppress a warning from LIME about feature discretization
warnings.filterwarnings("ignore", category=UserWarning, module='lime')

# Define the model function
lookup_table = {1.0: 1.0, 0.0: 0.0}
def f(inputs):
    """
    Model function that processes a batch of inputs for LIME.
    LIME passes a numpy array of shape (n_samples, n_features).
    """
    predictions = []
    for row in inputs:
        input1, input2 = row[0], row[1]
        # The default case depends only on input2
        default_value = input1 * 0 + input2 * 0.5 + 0.5
        # The get() method retrieves the value or returns the default
        pred = lookup_table.get(float(input1), default_value)
        predictions.append(pred)
    # LIME regression mode expects a 1D array of predictions
    return np.array(predictions)

# Define baseline data and feature names
# This data establishes that input1 values of 0 and 1 are "known" categories.
training_data = np.array([[0., 10.], [1., 20.]])
feature_names = ['input1', 'input2']

# --- Case i): E = (0.0, 0.0) ---
# For an "in-distribution" point, we simulate LIME's default behavior
# of treating features with few unique values in training as categorical.
print("--- Case i) E = (0.0, 0.0) ---")
explicand_1 = np.array([0.0, 0.0])

# Use discretize_continuous=True (LIME's default)
explainer_1 = lime.lime_tabular.LimeTabularExplainer(
    training_data=training_data,
    feature_names=feature_names,
    mode='regression',
    discretize_continuous=True
)

# Generate the explanation
exp1 = explainer_1.explain_instance(explicand_1, f, num_features=2)
weights1 = exp1.as_list()

print("Explanation for E = (0.0, 0.0):")
print("Feature Importances (Weights):")
importance1 = {}
for feature, weight in weights1:
    # The feature name is the first word in the description string
    feature_name = feature.split(' ')[0]
    importance1[feature_name] = weight
    print(f"  {feature}: {weight:.4f}")

# Compare absolute importances
abs_imp1_input1 = abs(importance1.get('input1', 0.0))
abs_imp1_input2 = abs(importance1.get('input2', 0.0))

print(f"\nAbsolute importance for input1: {abs_imp1_input1:.4f}")
print(f"Absolute importance for input2: {abs_imp1_input2:.4f}")
if abs_imp1_input1 > abs_imp1_input2:
    print("Result: For E=(0.0, 0.0), input1 is more important.")
else:
    print("Result: For E=(0.0, 0.0), input2 is more important.")

print("\n" + "="*40 + "\n")

# --- Case ii): E = (-1.0, -1.0) ---
# For an "out-of-distribution" point, the local continuous behavior of the model
# is more relevant. We simulate this by setting discretize_continuous=False.
print("--- Case ii) E = (-1.0, -1.0) ---")
explicand_2 = np.array([-1.0, -1.0])

# Use discretize_continuous=False to analyze local continuous behavior
explainer_2 = lime.lime_tabular.LimeTabularExplainer(
    training_data=training_data,
    feature_names=feature_names,
    mode='regression',
    discretize_continuous=False,
    sample_around_instance=True # focuses perturbations locally
)

# Generate the explanation
exp2 = explainer_2.explain_instance(explicand_2, f, num_features=2)
weights2 = exp2.as_list()

print("Explanation for E = (-1.0, -1.0):")
print("Feature Importances (Weights):")
importance2 = {}
for feature, weight in weights2:
    feature_name = feature.split(' ')[0]
    importance2[feature_name] = weight
    print(f"  {feature}: {weight:.4f}")

# Compare absolute importances
abs_imp2_input1 = abs(importance2.get('input1', 0.0))
abs_imp2_input2 = abs(importance2.get('input2', 0.0))

print(f"\nAbsolute importance for input1: {abs_imp2_input1:.4f}")
print(f"Absolute importance for input2: {abs_imp2_input2:.4f}")
if abs_imp2_input1 > abs_imp2_input2:
    print("Result: For E=(-1.0, -1.0), input1 is more important.")
else:
    print("Result: For E=(-1.0, -1.0), input2 is more important.")
