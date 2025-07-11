import numpy as np
import lime
import lime.lime_tabular
import warnings

# Suppress a known warning from LIME that is not relevant here
warnings.filterwarnings("ignore", category=UserWarning, module='lime')

def solve():
    """
    Generates LIME explanations for a given model and two input cases to determine feature importance.
    """
    # 1. Define the model function f for LIME
    lookup_table = {1.0: 1.0, 0.0: 0.0}

    def f_for_lime(numpy_array):
        """
        Prediction function that takes a numpy array and returns predictions.
        """
        n_samples = numpy_array.shape[0]
        predictions = np.zeros(n_samples)
        for i in range(n_samples):
            input1 = numpy_array[i, 0]
            input2 = numpy_array[i, 1]
            # The model logic provided in the problem
            predictions[i] = lookup_table.get(input1, input1 * 0 + input2 * 0.5 + 0.5)
        return predictions

    # 2. Define baseline data and explainer with default hyperparams
    # We interpret "baseline dataset is the same as the lookup table" as data
    # where features take on values 0 and 1.
    baseline_data = np.array([[0.0, 0.0], [1.0, 1.0], [0.0, 1.0], [1.0, 0.0]])
    feature_names = ['input1', 'input2']

    explainer = lime.lime_tabular.LimeTabularExplainer(
        training_data=baseline_data,
        feature_names=feature_names,
        mode='regression'
        # We use default hyperparams, including num_samples=5000 and discretize_continuous=True
    )

    # 3. Analyze Case (i): E is in the baseline dataset
    print("--- Case (i): Explicand E = (0.0, 0.0) ---")
    explicand_i = np.array([0.0, 0.0])
    exp_i = explainer.explain_instance(explicand_i, f_for_lime, num_features=2)
    
    # Construct and print the local model equation
    intercept_i = exp_i.intercept[0]
    rules_i = exp_i.as_list()
    equation_parts_i = [f"({weight:.4f} * I({rule}))" for rule, weight in rules_i]
    print(f"Local Model: prediction ≈ {' + '.join(equation_parts_i)} + {intercept_i:.4f}")

    # Determine most important feature for Case i
    if not rules_i or abs(rules_i[0][1]) > (abs(rules_i[1][1]) if len(rules_i) > 1 else 0):
        # as_list is sorted by importance, so check the first element
        most_important_i = 'input1' if 'input1' in rules_i[0][0] else 'input2'
    else:
        most_important_i = 'input2' if 'input2' in rules_i[1][0] else 'input1'
    print(f"Conclusion: {most_important_i} is more important for Case (i).\n")


    # 4. Analyze Case (ii): E is NOT in the baseline dataset
    print("--- Case (ii): Explicand E = (-1.0, -1.0) ---")
    explicand_ii = np.array([-1.0, -1.0])
    exp_ii = explainer.explain_instance(explicand_ii, f_for_lime, num_features=2)
    
    # Construct and print the local model equation
    intercept_ii = exp_ii.intercept[0]
    rules_ii = exp_ii.as_list()
    equation_parts_ii = [f"({weight:.4f} * I({rule}))" for rule, weight in rules_ii]
    print(f"Local Model: prediction ≈ {' + '.join(equation_parts_ii)} + {intercept_ii:.4f}")

    # Determine most important feature for Case ii
    if not rules_ii or abs(rules_ii[0][1]) > (abs(rules_ii[1][1]) if len(rules_ii) > 1 else 0):
        most_important_ii = 'input1' if 'input1' in rules_ii[0][0] else 'input2'
    else:
        most_important_ii = 'input2' if 'input2' in rules_ii[1][0] else 'input1'
    print(f"Conclusion: {most_important_ii} is more important for Case (ii).")


# Execute the analysis
solve()