import numpy as np
import lime
import lime.lime_tabular
import warnings

# Suppress potential warnings from lime regarding predict_proba
warnings.filterwarnings("ignore", category=UserWarning, module='lime')

def solve():
    """
    Solves the LIME feature importance problem.
    """
    # 1. Define the model function
    lookup_table = {1.0: 1.0, 0.0: 0.0}
    def f(data):
        """
        The model to be explained.
        It expects a numpy array of shape (n_samples, n_features).
        """
        predictions = np.zeros(data.shape[0])
        for i in range(data.shape[0]):
            input1 = data[i, 0]
            input2 = data[i, 1]
            if input1 in lookup_table:
                predictions[i] = lookup_table[input1]
            else:
                predictions[i] = input1 * 0 + input2 * 0.5 + 0.5
        return predictions

    # 2. Define the baseline dataset
    # Based on the lookup table keys {0.0, 1.0}
    baseline_data = np.array([[0.0, 0.0], [1.0, 1.0], [0.0, 1.0], [1.0, 0.0]], dtype=np.float64)
    feature_names = ['input1', 'input2']

    # 3. Set up the LIME explainer
    explainer = lime.lime_tabular.LimeTabularExplainer(
        training_data=baseline_data,
        feature_names=feature_names,
        class_names=['output'],
        mode='regression',
        verbose=False,
        discretize_continuous=False # Treat features as continuous numbers
    )

    # 4. Analyze Case i: E = (0.0, 0.0)
    print("--- Case i: Explaining E = (0.0, 0.0) ---")
    explicand_i = np.array([0.0, 0.0], dtype=np.float64)
    explanation_i = explainer.explain_instance(
        explicand_i,
        f,
        num_features=2,
        num_samples=5000
    )
    
    # Extract weights and intercept
    intercept_i = explanation_i.intercept['output']
    weights_i = dict(explanation_i.as_list())
    w1_i = weights_i.get('input1', 0.0)
    w2_i = weights_i.get('input2', 0.0)

    print("LIME Explanation as a local linear model:")
    print(f"Prediction ≈ {w1_i:.4f} * input1 + {w2_i:.4f} * input2 + {intercept_i:.4f}")

    if abs(w1_i) > abs(w2_i):
        important_feature_i = 'input1'
    else:
        important_feature_i = 'input2'
    print(f"Conclusion for Case i: The more important feature is '{important_feature_i}'.\n")


    # 5. Analyze Case ii: E = (-1.0, -1.0)
    print("--- Case ii: Explaining E = (-1.0, -1.0) ---")
    explicand_ii = np.array([-1.0, -1.0], dtype=np.float64)
    explanation_ii = explainer.explain_instance(
        explicand_ii,
        f,
        num_features=2,
        num_samples=5000
    )
    
    # Extract weights and intercept
    intercept_ii = explanation_ii.intercept['output']
    weights_ii = dict(explanation_ii.as_list())
    w1_ii = weights_ii.get('input1', 0.0)
    w2_ii = weights_ii.get('input2', 0.0)

    print("LIME Explanation as a local linear model:")
    print(f"Prediction ≈ {w1_ii:.4f} * input1 + {w2_ii:.4f} * input2 + {intercept_ii:.4f}")
    
    if abs(w1_ii) > abs(w2_ii):
        important_feature_ii = 'input1'
    else:
        important_feature_ii = 'input2'
    print(f"Conclusion for Case ii: The more important feature is '{important_feature_ii}'.")

    # Final result determination based on script output
    if important_feature_i == 'input1' and important_feature_ii == 'input2':
        final_answer = 'A'
    elif important_feature_i == 'input2' and important_feature_ii == 'input1':
        final_answer = 'B'
    elif important_feature_i == 'input1' and important_feature_ii == 'input1':
        final_answer = 'D'
    elif important_feature_i == 'input2' and important_feature_ii == 'input2':
        final_answer = 'E'
    else:
        final_answer = 'C' # Cannot say
        
    print(f"\nFinal conclusion matches answer choice {final_answer}.")
    print(f"<<<{final_answer}>>>")


if __name__ == '__main__':
    solve()