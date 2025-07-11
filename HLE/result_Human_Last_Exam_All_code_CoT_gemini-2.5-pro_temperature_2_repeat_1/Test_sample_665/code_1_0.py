import numpy as np
import lime
import lime.lime_tabular
import warnings

# Suppress a harmless warning from the LIME library regarding feature discretization
warnings.filterwarnings("ignore", category=UserWarning, module='lime')


def get_feature_importance(explanation):
    """Extracts feature names and their importance scores from a LIME explanation."""
    importance_map = {}
    for feature, weight in explanation.as_list():
        # Feature name is the first part, e.g., "input1" in "input1 <= 0.00"
        feature_name = feature.split(' ')[0]
        importance_map[feature_name] = abs(weight)
    return importance_map

# The model definition provided by the user
lookup_table = {1.0: 1.0, 0.0: 0.0}
def f(input_array):
    """
    The black-box model function.
    LIME's explainer expects a function that takes a numpy array (n_samples, n_features)
    and returns a numpy array (n_samples,).
    """
    predictions = []
    for row in input_array:
        input1, input2 = row
        # Use np.isclose for robust float comparisons
        if np.isclose(input1, 1.0):
            pred = 1.0
        elif np.isclose(input1, 0.0):
            pred = 0.0
        else:
            pred = input1 * 0 + input2 * 0.5 + 0.5
        predictions.append(pred)
    return np.array(predictions)

# Baseline dataset setup as per the problem description.
# "baseline dataset is the same as the lookup table" is interpreted as
# input1 values are drawn from the keys of the lookup table.
np.random.seed(42) # for reproducibility
train_input1 = np.random.choice(list(lookup_table.keys()), size=100)
train_input2 = np.random.rand(100) # Assuming a standard [0,1] distribution for the other feature
training_data = np.vstack([train_input1, train_input2]).T
feature_names = ['input1', 'input2']

# Initialize the LIME Explainer for tabular data
explainer = lime.lime_tabular.LimeTabularExplainer(
    training_data=training_data,
    feature_names=feature_names,
    class_names=['prediction'],
    mode='regression'
)

# --- Case i) E that belongs to the baseline dataset ---
explicand_i = np.array([0.0, 0.0])
explanation_i = explainer.explain_instance(
    explicand_i,
    f,
    num_features=2,
    num_samples=5000
)

print("--- Case (i): Explanation for E = (0.0, 0.0) ---")
importance_i = get_feature_importance(explanation_i)
print("Feature Importances:", importance_i)
if importance_i.get('input1', 0) > importance_i.get('input2', 0):
    print("Result: input1 is more important.")
else:
    print("Result: input2 is more important.")
print("\nLIME's raw explanation:")
for feature, weight in explanation_i.as_list():
    print(f"'{feature}' with weight {weight:.4f}")

print("\n" + "="*50 + "\n")

# --- Case ii) E that does not belong to the baseline dataset ---
explicand_ii = np.array([-1.0, -1.0])
explanation_ii = explainer.explain_instance(
    explicand_ii,
    f,
    num_features=2,
    num_samples=5000
)

print("--- Case (ii): Explanation for E = (-1.0, -1.0) ---")
importance_ii = get_feature_importance(explanation_ii)
print("Feature Importances:", importance_ii)
if importance_ii.get('input1', 0) > importance_ii.get('input2', 0):
    print("Result: input1 is more important.")
else:
    print("Result: input2 is more important.")
print("\nLIME's raw explanation:")
for feature, weight in explanation_ii.as_list():
    print(f"'{feature}' with weight {weight:.4f}")

<<<A>>>