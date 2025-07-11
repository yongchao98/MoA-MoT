import numpy as np
import lime
import lime.lime_tabular
import warnings

# Suppress a warning from LIME about feature discretization
warnings.filterwarnings("ignore", category=UserWarning, module='lime')

# Define the black-box model f
lookup_table = {1.0: 1.0, 0.0: 0.0}

def f(input_array):
    """
    The model function to be explained.
    LIME for classification expects (n_samples, n_classes) output.
    We'll treat the output as the probability of "Class 1".
    """
    input1_col = input_array[:, 0]
    input2_col = input_array[:, 1]
    
    # Initialize output array with the default case
    output = input1_col * 0 + input2_col * 0.5 + 0.5
    
    # Apply the lookup table logic using masks
    mask0 = (input1_col == 0.0)
    mask1 = (input1_col == 1.0)
    
    if np.any(mask0):
        output[mask0] = lookup_table[0.0]
    if np.any(mask1):
        output[mask1] = lookup_table[1.0]
    
    # LIME's classification mode needs probabilities for each class
    # We return [P(Class 0), P(Class 1)]
    return np.vstack([1 - output, output]).T

# The problem states "baseline dataset is the same as the lookup table".
# This means the training data for input1 primarily consists of 0.0 and 1.0.
# We will use a minimal dataset to represent this for the explainer's statistics.
training_data = np.array([[0.0, 0.0], [1.0, 1.0]])

# Create a LIME Tabular Explainer
explainer = lime.lime_tabular.LimeTabularExplainer(
    training_data=training_data,
    feature_names=['input1', 'input2'],
    class_names=['Class 0', 'Class 1'],
    mode='classification',
    random_state=42,
    discretize_continuous=False # Use raw features for clearer coefficient interpretation
)

# Case i): Explicand E belongs to the baseline dataset
e1 = np.array([0.0, 0.0])
exp1 = explainer.explain_instance(e1, f, num_features=2)
pred_e1 = f(e1.reshape(1, -1))[0, 1]

print("--- Case i): E = (0.0, 0.0) ---")
print(f"Model prediction P(Class 1) for E1: {pred_e1:.4f}")
print("Feature importances (coefficients of local linear model):")
imp1 = sorted(exp1.as_list(label=1), key=lambda x: abs(x[1]), reverse=True)
for feature, weight in imp1:
    print(f"- {feature}: {weight:.4f}")
print("Conclusion for i): 'input1' is more important due to the model's discontinuity at this point.")
print("\n" + "="*50 + "\n")


# Case ii): Explicand E does not belong to the baseline dataset
e2 = np.array([-1.0, -1.0])
exp2 = explainer.explain_instance(e2, f, num_features=2)
pred_e2 = f(e2.reshape(1, -1))[0, 1]

print("--- Case ii): E = (-1.0, -1.0) ---")
print(f"Model prediction P(Class 1) for E2: {pred_e2:.4f}")
print("Feature importances (coefficients of local linear model):")
imp2 = sorted(exp2.as_list(label=1), key=lambda x: abs(x[1]), reverse=True)
for feature, weight in imp2:
    print(f"- {feature}: {weight:.4f}")
print("Conclusion for ii): 'input2' is more important as the model is locally linear and only depends on input2.")
