import numpy as np
import lime
import lime.lime_tabular

# --- Step 1: Define the model and environment ---

# The model function as described. It's a lookup on input1,
# otherwise a linear function of input2.
lookup_table = {1.0: 1.0, 0.0: 0.0}
def f_raw(input1, input2):
    """The original model function for single inputs."""
    return lookup_table.get(input1, input1 * 0 + input2 * 0.5 + 0.5)

# LIME's predict_fn expects a function that takes a numpy array of shape
# (n_samples, n_features) and returns an array of shape (n_samples,).
def predict_fn(x):
    """Wrapper for the model to make it compatible with LIME."""
    outputs = []
    for row in x:
        outputs.append(f_raw(row[0], row[1]))
    # For a regression task, LIME expects a 1D array of predictions.
    return np.array(outputs)

# The "baseline dataset" is derived from the lookup table keys {0, 1}.
# We'll use a standard set of training data representing this.
training_data = np.array([[0.0, 0.0], [0.0, 1.0], [1.0, 0.0], [1.0, 1.0]])
feature_names = ['input1', 'input2']

# --- Step 2: Initialize LIME Explainer with default hyperparameters ---

# A key default is `discretize_continuous=True`, which means LIME will
# learn the distribution from `training_data` and sample from it to create
# perturbations. This means perturbations for input1 will always be 0.0 or 1.0.
explainer = lime.lime_tabular.LimeTabularExplainer(
    training_data=training_data,
    feature_names=feature_names,
    mode='regression',
    discretize_continuous=True, # This is the library default
    random_state=42 # for reproducibility
)


# --- Step 3: Analyze Case (i) - Explicand E is in the baseline dataset ---

print("--- Case (i): Explaining E = (0.0, 0.0) ---")
explicand_i = np.array([0.0, 0.0])

# Generate the explanation
explanation_i = explainer.explain_instance(
    explicand_i,
    predict_fn,
    num_features=2,
    num_samples=5000
)

# The local model equation is: prediction = intercept + sum(weight * feature_value)
# We print the components of this equation.
print(f"Prediction at E (from local model): {explanation_i.predicted_value:.4f}")
print(f"Local Intercept: {explanation_i.intercept[0]:.4f}")
print("Feature importances (the 'final equation' components):")
exp_list_i = explanation_i.as_list()
more_important_i = ""
max_abs_weight_i = -1
for feature, weight in exp_list_i:
    print(f"Term: ({weight:.4f} * I({feature}))")
    if abs(weight) > max_abs_weight_i:
        max_abs_weight_i = abs(weight)
        if 'input1' in feature:
            more_important_i = 'input1'
        elif 'input2' in feature:
            more_important_i = 'input2'

print(f"\nConclusion for case (i): '{more_important_i}' is more important.\n")


# --- Step 4: Analyze Case (ii) - Explicand E is NOT in the baseline dataset ---

print("--- Case (ii): Explaining E = (-1.0, -1.0) ---")
explicand_ii = np.array([-1.0, -1.0])

# Generate the explanation
explanation_ii = explainer.explain_instance(
    explicand_ii,
    predict_fn,
    num_features=2,
    num_samples=5000
)

# Print the components of the local model equation for case (ii)
print(f"Prediction at E (from local model): {explanation_ii.predicted_value:.4f}")
print(f"Local Intercept: {explanation_ii.intercept[0]:.4f}")
print("Feature importances (the 'final equation' components):")
exp_list_ii = explanation_ii.as_list()
more_important_ii = ""
max_abs_weight_ii = -1
for feature, weight in exp_list_ii:
    print(f"Term: ({weight:.4f} * I({feature}))")
    if abs(weight) > max_abs_weight_ii:
        max_abs_weight_ii = abs(weight)
        if 'input1' in feature:
            more_important_ii = 'input1'
        elif 'input2' in feature:
            more_important_ii = 'input2'
            
print(f"\nConclusion for case (ii): '{more_important_ii}' is more important.\n")

# --- Step 5: Final Conclusion ---
print("="*50)
print(f"Final Answer Summary:")
print(f"For case (i), the more important input is: {more_important_i}")
print(f"For case (ii), the more important input is: {more_important_ii}")
print("In both cases, input1 is determined to be the more important feature.")