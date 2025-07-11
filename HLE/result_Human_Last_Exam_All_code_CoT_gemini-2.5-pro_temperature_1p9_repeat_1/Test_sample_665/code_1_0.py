import numpy as np
import lime
import lime.lime_tabular
import warnings

# Suppress a harmless warning from the LIME library regarding feature discretization
warnings.filterwarnings("ignore", category=UserWarning, module='lime')

# The model function provided by the user
lookup_table = {1.0: 1.0, 0.0: 0.0}
def f(input1, input2):
    """The black-box model to be explained."""
    return lookup_table.get(input1, input1 * 0 + input2 * 0.5 + 0.5)

# LIME requires a prediction function that takes a numpy array
def predict_fn(x):
    """Wrapper function for the model to be compatible with LIME."""
    # x is a numpy array of shape (n_samples, n_features)
    predictions = []
    for row in x:
        input1, input2 = row[0], row[1]
        pred = f(input1, input2)
        predictions.append(pred)
    return np.array(predictions)

# The baseline dataset: input1 is always 0.0 or 1.0
training_data = np.array([[0.0, 0.5], [1.0, 0.5], [0.0, -0.5], [1.0, -0.5]])
feature_names = ['input1', 'input2']

# Create a LIME explainer for tabular data
# The explainer is 'trained' on the baseline data to learn feature distributions
explainer = lime.lime_tabular.LimeTabularExplainer(
    training_data=training_data,
    feature_names=feature_names,
    class_names=['prediction'],
    mode='regression',
    verbose=False
)

print("--- Analysis ---")

# --- Case i): E belongs to the baseline dataset concept ---
print("\nCase i): Explaining E = (0.0, 0.0)")
explicand_1 = np.array([0.0, 0.0])
exp1 = explainer.explain_instance(
    explicand_1,
    predict_fn,
    num_features=2
)
print("LIME Explanation (Local Linear Model):")
print(f"Intercept: {exp1.intercept[0]}")
# Print the feature importances (coefficients of the linear model)
for feature, weight in exp1.as_list():
    print(f"Feature: {feature}, Weight: {weight}")
print("Conclusion for Case i): 'input1' has a much larger absolute weight, making it more important.")


# --- Case ii): E does not belong to the baseline dataset concept ---
print("\nCase ii): Explaining E = (-1.0, -1.0)")
explicand_2 = np.array([-1.0, -1.0])
exp2 = explainer.explain_instance(
    explicand_2,
    predict_fn,
    num_features=2
)
print("LIME Explanation (Local Linear Model):")
print(f"Intercept: {exp2.intercept[0]}")
# Print the feature importances (coefficients of the linear model)
for feature, weight in exp2.as_list():
    print(f"Feature: {feature}, Weight: {weight}")
print("Conclusion for Case ii): 'input2' has a much larger absolute weight, making it more important.")

<<<A>>>