import numpy as np

# Define the model function
lookup_table = {1.0: 1.0, 0.0: 0.0}
def f(input1, input2):
    """The black-box model to be explained."""
    # Convert inputs to float, as they might be numpy floats
    input1 = float(input1)
    if input1 in lookup_table:
        return lookup_table[input1]
    else:
        return input1 * 0 + input2 * 0.5 + 0.5

def explain_local_behavior(explicand, model, delta=1e-6):
    """
    Analyzes the local behavior of a model around an explicand
    by observing the change in output for small perturbations in each input.
    This simulates the core idea of a local explanation like LIME.
    """
    E1, E2 = explicand
    base_prediction = model(E1, E2)

    # Perturb input1 and calculate change
    pred_on_pert1 = model(E1 + delta, E2)
    importance1 = abs(pred_on_pert1 - base_prediction)

    # Perturb input2 and calculate change
    pred_on_pert2 = model(E1, E2 + delta)
    importance2 = abs(pred_on_pert2 - base_prediction)
    
    return importance1, importance2

# --- Setup ---
delta = 1e-6
# Case i) E belongs to the baseline dataset, e.g., (0.0, 0.0)
E_i = (0.0, 0.0)
imp1_i, imp2_i = explain_local_behavior(E_i, f, delta)

# Case ii) E does not belong to the baseline dataset, e.g., (-1.0, -1.0)
E_ii = (-1.0, -1.0)
imp1_ii, imp2_ii = explain_local_behavior(E_ii, f, delta)

print("Analysis of local feature importance:")
print("-" * 50)

# --- Analysis for Case i ---
print("Case i: Explicand E = (0.0, 0.0)")
base_pred_i = f(E_i[0], E_i[1])
pert1_pred_i = f(E_i[0] + delta, E_i[1])
pert2_pred_i = f(E_i[0], E_i[1] + delta)

print(f"Prediction for E: f({E_i[0]}, {E_i[1]}) = {base_pred_i}")
print(f"Equation for input1 importance: |f({E_i[0] + delta}, {E_i[1]}) - f({E_i[0]}, {E_i[1]})| = |{pert1_pred_i:.6f} - {base_pred_i}| = {imp1_i:.6f}")
print(f"Equation for input2 importance: |f({E_i[0]}, {E_i[1] + delta}) - f({E_i[0]}, {E_i[1]})| = |{pert2_pred_i} - {base_pred_i}| = {imp2_i:.6f}")

if imp1_i > imp2_i:
    result_i = "input1 is more important"
else:
    result_i = "input2 is more important"
print(f"Conclusion for Case i: {result_i}")
print("-" * 50)

# --- Analysis for Case ii ---
print("Case ii: Explicand E = (-1.0, -1.0)")
base_pred_ii = f(E_ii[0], E_ii[1])
pert1_pred_ii = f(E_ii[0] + delta, E_ii[1])
pert2_pred_ii = f(E_ii[0], E_ii[1] + delta)

print(f"Prediction for E: f({E_ii[0]}, {E_ii[1]}) = {base_pred_ii}")
print(f"Equation for input1 importance: |f({E_ii[0] + delta}, {E_ii[1]}) - f({E_ii[0]}, {E_ii[1]})| = |{pert1_pred_ii:.6f} - {base_pred_ii}| = {imp1_ii:.6f}")
print(f"Equation for input2 importance: |f({E_ii[0]}, {E_ii[1] + delta}) - f({E_ii[0]}, {E_ii[1]})| = |{pert2_pred_ii:.6f} - {base_pred_ii}| = {imp2_ii:.6f}")

if imp1_ii > imp2_ii:
    result_ii = "input1 is more important"
else:
    result_ii = "input2 is more important"
print(f"Conclusion for Case ii: {result_ii}")
print("-" * 50)