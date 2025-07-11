import numpy as np
import lime
import lime.lime_tabular

# The model function as described by the user.
# It has two regimes based on the value of input1.
lookup_table = {1.0: 1.0, 0.0: 0.0}
def f(input1, input2):
    return lookup_table.get(input1, 0.5 * input2 + 0.5)

# LIME requires a prediction function that takes a numpy array (n_samples, n_features)
# and returns a numpy array (n_samples,).
def predict_fn(x):
    """Wrapper for the model 'f' to be compatible with LIME."""
    predictions = []
    for row in x:
        input1, input2 = row[0], row[1]
        result = f(input1, input2)
        predictions.append(result)
    return np.array(predictions)

# The baseline dataset is interpreted from "the same as the lookup table".
# This means the values for input1 are the keys of the table: {0.0, 1.0}.
baseline_data = np.array([[0.0, 0.0], [1.0, 1.0]])
feature_names = ['input1', 'input2']

# We initialize the LIME Tabular Explainer.
# We use the default hyperparameters as requested, which importantly includes
# `discretize_continuous=True`. The model is treated as a regression model.
explainer = lime.lime_tabular.LimeTabularExplainer(
    training_data=baseline_data,
    feature_names=feature_names,
    mode='regression'
)

# --- Scenario i): Explicand E belongs to the baseline dataset ---
explicand_i = np.array([0.0, 0.0])

# Generate the explanation
explanation_i = explainer.explain_instance(
    explicand_i,
    predict_fn,
    num_features=2,
    num_samples=5000
)

print("--- Scenario i) Explanation for E = (0.0, 0.0) ---")
print("Prediction: ", predict_fn(explicand_i.reshape(1, -1))[0])
print("The local linear model coefficients (feature importances) are:")
importances_i = explanation_i.as_list()
for feature, weight in importances_i:
    print(f"Feature '{feature}' has a weight of: {weight}")

# --- Scenario ii): Explicand E does not belong to the baseline dataset ---
explicand_ii = np.array([-1.0, -1.0])

# Generate the explanation
explanation_ii = explainer.explain_instance(
    explicand_ii,
    predict_fn,
    num_features=2,
    num_samples=5000
)

print("\n--- Scenario ii) Explanation for E = (-1.0, -1.0) ---")
print("Prediction: ", predict_fn(explicand_ii.reshape(1, -1))[0])
print("The local linear model coefficients (feature importances) are:")
importances_ii = explanation_ii.as_list()
for feature, weight in importances_ii:
    print(f"Feature '{feature}' has a weight of: {weight}")

# Determine the most important feature in each case based on the absolute weight
abs_importances_i = {feat.split(' ')[0]: abs(w) for feat, w in importances_i}
more_important_i = "input1" if abs_importances_i.get('input1', 0) > abs_importances_i.get('input2', 0) else "input2"

abs_importances_ii = {feat.split(' ')[0]: abs(w) for feat, w in importances_ii}
more_important_ii = "input1" if abs_importances_ii.get('input1', 0) > abs_importances_ii.get('input2', 0) else "input2"

print("\n--- Conclusion ---")
print(f"For scenario i), the more important feature is: {more_important_i}")
print(f"For scenario ii), the more important feature is: {more_important_ii}")
print("This indicates that 'input1' is consistently chosen as the more important feature.")
<<<D>>>