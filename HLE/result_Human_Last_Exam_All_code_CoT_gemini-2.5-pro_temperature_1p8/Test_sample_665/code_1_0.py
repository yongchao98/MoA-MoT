import lime
import lime.lime_tabular
import numpy as np
import warnings

# Suppress a harmless future warning from LIME
warnings.filterwarnings("ignore", category=FutureWarning)

def solve():
    """
    This function sets up the LIME explanation scenarios as described in the problem
    and prints the feature importances to determine the correct answer.
    """
    # Step 1: Define the model function
    # LIME's tabular explainer requires a function that takes a numpy array
    # of shape (n_samples, n_features) and returns a numpy array of shape (n_samples,).
    lookup_table = {1.0: 1.0, 0.0: 0.0}
    def f_model(x):
        """
        Vectorized model function for LIME.
        - If input1 is 1.0, output is 1.0.
        - If input1 is 0.0, output is 0.0.
        - Otherwise, output is 0.5 * input2 + 0.5.
        """
        outputs = np.zeros(x.shape[0])
        for i, (input1, input2) in enumerate(x):
            if input1 == 1.0:
                outputs[i] = 1.0
            elif input1 == 0.0:
                outputs[i] = 0.0
            else:
                outputs[i] = 0.5 * input2 + 0.5
        return outputs

    # Step 2: Define the baseline dataset
    # "baseline dataset is the same as the lookup table" is interpreted as the
    # set of points formed by the keys/values. A robust training set from
    # these values is used for LIME's perturbation sampling.
    baseline_data = np.array([
        [0.0, 0.0],
        [1.0, 1.0],
        [0.0, 1.0],
        [1.0, 0.0]
    ])
    feature_names = ['input1', 'input2']

    # Step 3: Create the LIME Tabular Explainer
    explainer = lime.lime_tabular.LimeTabularExplainer(
        training_data=baseline_data,
        feature_names=feature_names,
        mode='regression',
        random_state=42 # for reproducible perturbations
    )

    # --- Scenario i) E belongs to the baseline dataset: E = (0.0, 0.0) ---
    explicand_i = np.array([0.0, 0.0])
    explanation_i = explainer.explain_instance(
        explicand_i,
        f_model,
        num_features=2
    )

    # Get the raw feature weights from the explanation
    exp_map_i = explanation_i.as_map()[1]
    weights_i = {feature_names[i]: w for i, w in exp_map_i}

    print("--- Scenario i: E = (0.0, 0.0) ---")
    print(f"Model prediction f(0.0, 0.0): {f_model(explicand_i.reshape(1,-1))[0]}")
    print("The local linear model LIME found is:")
    print(f"Prediction ≈ {weights_i.get('input1', 0):.4f} * [input1] + {weights_i.get('input2', 0):.4f} * [input2] + intercept")
    if abs(weights_i.get('input1', 0)) > abs(weights_i.get('input2', 0)):
        print("Conclusion: input1 is more important for case i).")
    else:
        print("Conclusion: input2 is more important for case i).")
    print("-" * 35)
    print()

    # --- Scenario ii) E does not belong to the baseline dataset: E = (-1.0, -1.0) ---
    explicand_ii = np.array([-1.0, -1.0])
    explanation_ii = explainer.explain_instance(
        explicand_ii,
        f_model,
        num_features=2
    )

    # Get the raw feature weights
    exp_map_ii = explanation_ii.as_map()[1]
    weights_ii = {feature_names[i]: w for i, w in exp_map_ii}

    print("--- Scenario ii: E = (-1.0, -1.0) ---")
    print(f"Model prediction f(-1.0, -1.0): {f_model(explicand_ii.reshape(1,-1))[0]}")
    print("The local linear model LIME found is:")
    print(f"Prediction ≈ {weights_ii.get('input1', 0):.4f} * [input1] + {weights_ii.get('input2', 0):.4f} * [input2] + intercept")
    if abs(weights_ii.get('input1', 0)) > abs(weights_ii.get('input2', 0)):
        print("Conclusion: input1 is more important for case ii).")
    else:
        print("Conclusion: input2 is more important for case ii).")
    print("-" * 35)

solve()