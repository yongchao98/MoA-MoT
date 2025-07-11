import numpy as np

# --- 1. Define placeholder values for each term in the model ---
# These values are chosen for illustrative purposes to show the structure of the equation.

# Learning and Regularization Terms
learning_rate_term = 0.01
mission_based_utility_term = 5.0
gradient_of_loss = 2.0
weight_regularization_term = 0.1
learning_utility_term = 0.5
decay_utility_term = 0.2
external_stimulus_impact_term = 1.0

# Pruning and Plasticity Terms
pruning_probability_term = 0.001
# For: Activation Function(-Utility-Based Pruning Term + Randomness Term)
utility_based_pruning_term = 3.0
pruning_randomness_term = 0.5
# For: Activation Function(|Weights|)
weights_magnitude = 4.0

# Global and Spatial Dynamics
global_randomness_term = 0.005
randomness_coefficient = 0.8
spatial_diffusion_term = 0.02

# Dynamic Threshold Terms (Negative contribution)
base_threshold = 1.0
fatigue_coefficient = 0.2
integral_recent_activity = 1.5 # ∫ from t-Δt to t [Recent Activity] dτ
cumulative_activity_coefficient = 0.05
integral_cumulative_activity = 10.0 # ∫ from 0 to t [Cumulative Activity] dτ

# Long-term Memory and Input Modulation
integral_memory_influence = 0.3 # ∫ from 0 to t [Memory Decay Term × Historical Influence] dτ
input_relevance_dropout = 0.1

# For simplicity, we assume Activation Function(x) = x
# This is a linear activation, in reality it would be a non-linear function like sigmoid or ReLU.

# --- 2. Calculate each part of the equation for Model A ---

term1 = learning_rate_term * (mission_based_utility_term + gradient_of_loss)

term2 = -learning_rate_term * (gradient_of_loss + weight_regularization_term)

term3 = -learning_rate_term * learning_utility_term * (gradient_of_loss + weight_regularization_term + decay_utility_term + external_stimulus_impact_term)

# Assuming linear activation for simplicity
term4 = -pruning_probability_term * (-utility_based_pruning_term + pruning_randomness_term)

term5 = -pruning_probability_term * weights_magnitude

term6 = global_randomness_term * randomness_coefficient

term7 = spatial_diffusion_term

# Note the internal signs: Fatigue (+) increases the subtracted threshold, Cumulative Activity (-) decreases it
term8 = -(base_threshold + fatigue_coefficient * integral_recent_activity - cumulative_activity_coefficient * integral_cumulative_activity)

term9 = integral_memory_influence

term10 = input_relevance_dropout

# --- 3. Calculate the final differential update ---
dw_dt = term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8 + term9 + term10

# --- 4. Print the full equation with the substituted numbers ---
print("Model A is the optimal choice. Its equation combines all the key features of a comprehensive neuromorphic model.")
print("Here is the calculation for ∂w(x, t) / ∂t based on the sample values:\n")
print(f"  {learning_rate_term} * ({mission_based_utility_term} + {gradient_of_loss})")
print(f"- {learning_rate_term} * ({gradient_of_loss} + {weight_regularization_term})")
print(f"- {learning_rate_term} * {learning_utility_term} * ({gradient_of_loss} + {weight_regularization_term} + {decay_utility_term} + {external_stimulus_impact_term})")
print(f"- {pruning_probability_term} * (-{utility_based_pruning_term} + {pruning_randomness_term})")
print(f"- {pruning_probability_term} * {weights_magnitude}")
print(f"+ {global_randomness_term} * {randomness_coefficient}")
print(f"+ {spatial_diffusion_term}")
print(f"- ({base_threshold} + {fatigue_coefficient} * {integral_recent_activity} - {cumulative_activity_coefficient} * {integral_cumulative_activity})")
print(f"+ {integral_memory_influence}")
print(f"+ {input_relevance_dropout}")

print("\n= " + " + ".join(f"({val:.5f})" for val in [term1, term2, term3, term4, term5, term6, term7, term8, term9, term10]))
print(f"\nFinal calculated value for ∂w(x, t) / ∂t = {dw_dt:.5f}")
