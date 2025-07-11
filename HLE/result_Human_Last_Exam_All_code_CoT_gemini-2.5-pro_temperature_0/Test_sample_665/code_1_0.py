import numpy as np
import lime
import lime.lime_tabular

# Step 1: Define the model and a wrapper for LIME
lookup_table = {1.0: 1.0, 0.0: 0.0}
def f(input1, input2):
    """The black-box model function."""
    # If input1 is a key in the lookup_table, the output depends only on input1.
    # Otherwise, the output depends only on input2.
    return lookup_table.get(input1, input1 * 0 + input2 * 0.5 + 0.5)

def predict_fn(x):
    """Wrapper function for LIME, which processes a batch of samples."""
    predictions = []
    for row in x:
        input1, input2 = row[0], row[1]
        predictions.append(f(input1, input2))
    return np.array(predictions)

# Step 2: Set up the LIME explainer
# The problem states "the baseline dataset is the same as the lookup table".
# We create a baseline dataset where input1 values are the keys of the lookup table.
baseline_dataset = np.array([[0.0, 0.0], [1.0, 1.0]])
feature_names = ['input1', 'input2']

# Create the explainer using default hyperparameters.
explainer = lime.lime_tabular.LimeTabularExplainer(
    training_data=baseline_dataset,
    feature_names=feature_names,
    mode='regression'
)

# Step 3: Analyze Case i, E that belongs to the baseline dataset
print("--- Analyzing Case i: E = (0.0, 0.0) ---")
instance_i = np.array([0.0, 0.0])
prediction_i = predict_fn(instance_i.reshape(1, -1))[0]
print(f"Model prediction for E=(0.0, 0.0): f({instance_i[0]}, {instance_i[1]}) = {prediction_i}")

explanation_i = explainer.explain_instance(
    instance_i,
    predict_fn,
    num_features=2
)
importances_i = explanation_i.as_list()

print("\nLIME Explanation for E=(0.0, 0.0):")
abs_importance_i = {'input1': 0, 'input2': 0}
for feature, weight in importances_i:
    if 'input1' in feature:
        abs_importance_i['input1'] += abs(weight)
    elif 'input2' in feature:
        abs_importance_i['input2'] += abs(weight)
    print(f"Feature: {feature}, Weight: {weight:.4f}")

most_important_i = 'input1' if abs_importance_i['input1'] > abs_importance_i['input2'] else 'input2'
print(f"\nConclusion for Case i: '{most_important_i}' is more important.")


# Step 4: Analyze Case ii, E that does not belong to the baseline dataset
print("\n--- Analyzing Case ii: E = (-1.0, -1.0) ---")
instance_ii = np.array([-1.0, -1.0])
prediction_ii = predict_fn(instance_ii.reshape(1, -1))[0]
print(f"Model prediction for E=(-1.0, -1.0): f({instance_ii[0]}, {instance_ii[1]}) = {prediction_ii}")

explanation_ii = explainer.explain_instance(
    instance_ii,
    predict_fn,
    num_features=2
)
importances_ii = explanation_ii.as_list()

print("\nLIME Explanation for E=(-1.0, -1.0):")
abs_importance_ii = {'input1': 0, 'input2': 0}
for feature, weight in importances_ii:
    if 'input1' in feature:
        abs_importance_ii['input1'] += abs(weight)
    elif 'input2' in feature:
        abs_importance_ii['input2'] += abs(weight)
    print(f"Feature: {feature}, Weight: {weight:.4f}")

most_important_ii = 'input1' if abs_importance_ii['input1'] > abs_importance_ii['input2'] else 'input2'
print(f"\nConclusion for Case ii: '{most_important_ii}' is more important.")

# Step 5: Final Answer
print("\n" + "="*40)
print(f"Summary:")
print(f"For Case i (E in baseline, e.g., (0.0, 0.0)), the most important feature is: {most_important_i}")
print(f"For Case ii (E not in baseline, e.g., (-1.0, -1.0)), the most important feature is: {most_important_ii}")
print("="*40)
<<<A>>>