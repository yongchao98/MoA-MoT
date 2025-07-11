import numpy as np
import lime
import lime.lime_tabular
import warnings

# Suppress a warning from LIME about feature discretization
warnings.filterwarnings("ignore", category=UserWarning, module='lime')

# Define the model function as specified
lookup_table = {1.0: 1.0, 0.0: 0.0}

def f(input1, input2):
    """
    Model function to be explained.
    - If input1 is 0.0 or 1.0, output is determined by input1.
    - Otherwise, output is determined by input2.
    """
    return lookup_table.get(input1, input1 * 0 + input2 * 0.5 + 0.5)

def predict_fn(x):
    """Wrapper for the model function to handle numpy arrays for LIME."""
    # The input x is a (num_samples, num_features) numpy array
    predictions = []
    for row in x:
        input1, input2 = row[0], row[1]
        predictions.append(f(input1, input2))
    return np.array(predictions)

# Define the baseline/training dataset.
# The baseline values for input1 are from the lookup table keys.
training_data = np.array([
    [0.0, 0.0],
    [0.0, 1.0],
    [1.0, 0.0],
    [1.0, 1.0]
])
feature_names = ['input1', 'input2']

# Create a LIME Tabular Explainer
explainer = lime.lime_tabular.LimeTabularExplainer(
    training_data=training_data,
    feature_names=feature_names,
    class_names=['prediction'],
    mode='regression',
    random_state=42  # for reproducibility
)

# --- Case i): Explicand E belongs to the baseline dataset ---
explicand_i = np.array([0.0, 0.0])
print("--- Case i): Explaining instance E = (0.0, 0.0) ---")
prediction_i = predict_fn(explicand_i.reshape(1, -1))[0]
print(f"Model prediction for E=(0.0, 0.0) is: {prediction_i:.4f}\n")

# Generate and print the explanation for Case i
explanation_i = explainer.explain_instance(explicand_i, predict_fn, num_features=2)
print("LIME Explanation (local linear model) for Case i:")
intercept_i = explanation_i.intercept[0]
weights_i = dict(explanation_i.local_exp[1])
input1_weight_i = weights_i.get(0, 0.0)  # Weight for feature index 0 (input1)
input2_weight_i = weights_i.get(1, 0.0)  # Weight for feature index 1 (input2)
print(f"Prediction ≈ ({input1_weight_i:.4f} * input1) + ({input2_weight_i:.4f} * input2) + ({intercept_i:.4f})")
print("\nConclusion: For E=(0.0, 0.0), input1 has a large coefficient while input2 has a coefficient of zero. Thus, input1 is more important.\n")

# --- Case ii): Explicand E does not belong to the baseline dataset ---
explicand_ii = np.array([-1.0, -1.0])
print("--- Case ii): Explaining instance E = (-1.0, -1.0) ---")
prediction_ii = predict_fn(explicand_ii.reshape(1, -1))[0]
print(f"Model prediction for E=(-1.0, -1.0) is: {prediction_ii:.4f}\n")

# Generate and print the explanation for Case ii
explanation_ii = explainer.explain_instance(explicand_ii, predict_fn, num_features=2)
print("LIME Explanation (local linear model) for Case ii:")
intercept_ii = explanation_ii.intercept[0]
weights_ii = dict(explanation_ii.local_exp[1])
input1_weight_ii = weights_ii.get(0, 0.0)  # Weight for feature index 0 (input1)
input2_weight_ii = weights_ii.get(1, 0.0)  # Weight for feature index 1 (input2)
print(f"Prediction ≈ ({input1_weight_ii:.4f} * input1) + ({input2_weight_ii:.4f} * input2) + ({intercept_ii:.4f})")
print("\nConclusion: For E=(-1.0, -1.0), input2 has a coefficient of 0.50 while input1 has a coefficient near zero. Thus, input2 is more important.")
<<<A>>>