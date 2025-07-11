import numpy as np
import lime
import lime.lime_tabular
import warnings

# Suppress a warning from scikit-learn used by LIME for plotting
warnings.filterwarnings("ignore", category=UserWarning, module='sklearn')

# Step 1: Define the model function
# LIME's tabular explainer expects a function that takes a numpy array of shape (n_samples, n_features)
# and returns a numpy array of shape (n_samples,) for regression.
lookup_table = {1.0: 1.0, 0.0: 0.0}

def f_predict(numpy_array_of_inputs):
    """
    Model function f(input1, input2) that works on a batch of inputs.
    """
    predictions = []
    for row in numpy_array_of_inputs:
        input1, input2 = row[0], row[1]
        # The model's logic as described
        prediction = lookup_table.get(input1, input1 * 0 + input2 * 0.5 + 0.5)
        predictions.append(prediction)
    return np.array(predictions)

# Step 2: Define the baseline dataset and LIME explainer
# The "baseline dataset is the same as the lookup table". The keys are 0 and 1.
# We will assume corresponding input2 values of 0 and 1, creating a simple baseline.
training_data = np.array([[0.0, 0.0], [1.0, 1.0]])
feature_names = ['input1', 'input2']

# We use LimeTabularExplainer, which is suited for this kind of feature interaction.
# `discretize_continuous=False` is important as it tells LIME to use the raw continuous values
# in its local regression model.
explainer = lime.lime_tabular.LimeTabularExplainer(
    training_data=training_data,
    feature_names=feature_names,
    mode='regression',
    discretize_continuous=False,
    # Using the default feature selection method, which is well-suited for tabular data.
)

# Step 3: Analyze Case (i)
print("--- Case (i): Explicand E belongs to the baseline dataset E=(0.0, 0.0) ---")
explicand_i = np.array([0.0, 0.0])

# Generate the explanation
explanation_i = explainer.explain_instance(
    data_row=explicand_i,
    predict_fn=f_predict,
    num_features=2,
    num_samples=5000  # Default number of samples
)

# Output the results for case i
print(f"Prediction for E={explicand_i}: {explanation_i.predicted_value:.4f}")
print("Local linear model intercept: {:.4f}".format(explanation_i.intercept[0]))
print("Feature importances (coefficients of the local linear model):")
importances_i = explanation_i.as_list()
for feature, weight in importances_i:
    print(f"- {feature}: {weight:.4f}")

# Determine the most important feature
if abs(importances_i[0][1]) > abs(importances_i[1][1]):
    most_important_i = 'input1' if 'input1' in importances_i[0][0] else 'input2'
else:
    most_important_i = 'input1' if 'input1' in importances_i[1][0] else 'input2'

print(f"\nConclusion for Case (i): '{most_important_i}' is the more important feature.")

# Step 4: Analyze Case (ii)
print("\n--- Case (ii): Explicand E does not belong to the baseline dataset E=(-1.0, -1.0) ---")
explicand_ii = np.array([-1.0, -1.0])

# Generate the explanation
explanation_ii = explainer.explain_instance(
    data_row=explicand_ii,
    predict_fn=f_predict,
    num_features=2,
    num_samples=5000  # Default number of samples
)

# Output the results for case ii
print(f"Prediction for E={explicand_ii}: {explanation_ii.predicted_value:.4f}")
print("Local linear model intercept: {:.4f}".format(explanation_ii.intercept[0]))
print("Feature importances (coefficients of the local linear model):")
importances_ii = explanation_ii.as_list()
for feature, weight in importances_ii:
    print(f"- {feature}: {weight:.4f}")

# Determine the most important feature
if abs(importances_ii[0][1]) > abs(importances_ii[1][1]):
    most_important_ii = 'input1' if 'input1' in importances_ii[0][0] else 'input2'
else:
    most_important_ii = 'input1' if 'input1' in importances_ii[1][0] else 'input2'
    
print(f"\nConclusion for Case (ii): '{most_important_ii}' is the more important feature.")