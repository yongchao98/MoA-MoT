# Plan:
# 1. Define plausible example variables for each term in Model A.
# 2. Construct the full equation using these variables.
# 3. Calculate the final result for the rate of change of the weight (dw_dt).
# 4. Print the full equation with the numerical values substituted.
# 5. Print the final calculated value.

# --- 1. Define variables for Model A ---
# These are example values to illustrate the calculation.
# Learning and Loss related terms
learning_rate_term = 0.01
mission_based_utility_term = 0.5
grad_loss_wrt_weights = 1.2
weight_regularization_term = 0.02

# Utility and stimulus related terms
learning_utility_term = 0.8
decay_utility_term = 0.005
external_stimulus_impact_term = 0.1

# Pruning related terms
pruning_probability_term = 0.001
# Assume activation function outputs some value. We use two examples.
act_func_val_1 = 0.9 # For Activation Function(− Utility-Based Pruning Term + Randomness Term)
act_func_val_2 = 0.7 # For Activation Function(|Weights|)

utility_based_pruning_term = 0.3
randomness_term = 0.05
abs_weights = 1.5

# Global effects
global_randomness_term = 0.0001
randomness_coefficient = -0.4  # Can be positive or negative
spatial_diffusion_term = 0.03

# Adaptive Threshold related terms (integrals are represented as pre-calculated values)
base_threshold = 0.2
fatigue_coefficient = 0.1
integral_recent_activity = 0.6
cumulative_activity_coefficient = 0.01
integral_cumulative_activity = 25.0

# Memory and Input related terms
integral_memory_influence = 0.08 # ∫ from 0 to t [Memory Decay Term × Historical Influence] dτ
input_relevance_term = 0.95
dropout_mask = 1.0 # 1 if connection is kept, 0 if dropped out

# --- 2. Calculate each part of the equation ---

# Term 1: Mission-driven update
term1 = learning_rate_term * (mission_based_utility_term + grad_loss_wrt_weights)

# Term 2: Regularized gradient descent (subtraction)
term2 = -learning_rate_term * (grad_loss_wrt_weights + weight_regularization_term)

# Term 3: Utility-driven learning modulation (subtraction)
term3 = -learning_rate_term * learning_utility_term * (grad_loss_wrt_weights + weight_regularization_term + decay_utility_term + external_stimulus_impact_term)

# Term 4: Utility-based pruning (subtraction)
term4 = -pruning_probability_term * act_func_val_1 # Where act_func_val_1 = Activation_Function(-utility_based_pruning_term + randomness_term)

# Term 5: Magnitude-based pruning (subtraction)
term5 = -pruning_probability_term * act_func_val_2 # Where act_func_val_2 = Activation_Function(|Weights|)

# Term 6: Global randomness (addition)
term6 = global_randomness_term * randomness_coefficient

# Term 7: Spatial diffusion (addition)
term7 = spatial_diffusion_term

# Term 8: Adaptive threshold (subtraction)
adaptive_threshold = (base_threshold + fatigue_coefficient * integral_recent_activity - cumulative_activity_coefficient * integral_cumulative_activity)
term8 = -adaptive_threshold

# Term 9: Memory influence (addition)
term9 = integral_memory_influence

# Term 10: Input relevance / dropout (addition)
term10 = input_relevance_term * dropout_mask

# --- 3. Calculate the final result ---
dw_dt = term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8 + term9 + term10


# --- 4. Print the full equation with numbers ---
print("Optimal Model A: Differential Updates (∂w/∂t) = ")
print(f"  {learning_rate_term} * ({mission_based_utility_term} + {grad_loss_wrt_weights})")
print(f"- {learning_rate_term} * ({grad_loss_wrt_weights} + {weight_regularization_term})")
print(f"- {learning_rate_term} * {learning_utility_term} * ({grad_loss_wrt_weights} + {weight_regularization_term} + {decay_utility_term} + {external_stimulus_impact_term})")
print(f"- {pruning_probability_term} * Activation_Function(-{utility_based_pruning_term} + {randomness_term})  # Evaluates to {act_func_val_1}")
print(f"- {pruning_probability_term} * Activation_Function(|{abs_weights}|) # Evaluates to {act_func_val_2}")
print(f"+ {global_randomness_term} * {randomness_coefficient}")
print(f"+ {spatial_diffusion_term}")
print(f"- ({base_threshold} + {fatigue_coefficient} * {integral_recent_activity} - {cumulative_activity_coefficient} * {integral_cumulative_activity})")
print(f"+ {integral_memory_influence}")
print(f"+ {input_relevance_term} * {dropout_mask}\n")

# --- 5. Print the final calculated value ---
print("--- Calculation ---")
print(f"∂w/∂t = {term1:.4f}  # Mission-Based Utility Term")
print(f"      + {term2:.4f}  # Regularized Gradient Term")
print(f"      + {term3:.4f}  # Learning Utility Term")
print(f"      + {term4:.4f}  # Utility-Based Pruning Term")
print(f"      + {term5:.4f}  # Magnitude-Based Pruning Term")
print(f"      + {term6:.4f}  # Global Randomness Term")
print(f"      + {term7:.4f}  # Spatial Diffusion Term")
print(f"      + {term8:.4f}  # Adaptive Threshold Term")
print(f"      + {term9:.4f}  # Memory Influence Term")
print(f"      + {term10:.4f}  # Input Relevance Term")
print("---------------------")
print(f"Final ∂w(x, t)/∂t = {dw_dt}")
<<<A>>>