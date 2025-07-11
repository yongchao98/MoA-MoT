import numpy as np

# --- Assign representative numerical values to each term in Model A ---

# Learning terms
learning_rate_term = 0.01
mission_based_utility_term = 5.0
gradient_of_loss = 3.0
weight_regularization_term = 0.1
learning_utility_term = 0.5
decay_utility_term = 0.2
external_stimulus_impact_term = 1.0

# Pruning terms (assuming linear activation for simplicity)
pruning_probability_term = 0.05
utility_based_pruning_term = 4.0
pruning_randomness_term = 0.5
weights_magnitude = 2.0

# Stochastic and spatial terms
global_randomness_term = 0.02
randomness_coefficient = 1.0
spatial_diffusion_term = 0.3

# Dynamic threshold and activity-dependent terms (integrals are simplified to representative values)
base_threshold = 0.2
fatigue_coefficient = 0.1
recent_activity_integral = 15.0
cumulative_activity_coefficient = 0.01
cumulative_activity_integral = 100.0

# Memory and input relevance terms
memory_decay_historical_influence_integral = 1.8 # ∫ [Memory Decay Term × Historical Influence] dτ
input_relevance_term = 0.7
dropout_mask = 1.0 # Neuron is active

# --- Calculate each component of the equation for Model A ---

# Term 1: Mission-driven learning
term1 = learning_rate_term * (mission_based_utility_term + gradient_of_loss)

# Term 2: Standard regularization
term2 = learning_rate_term * (gradient_of_loss + weight_regularization_term)

# Term 3: Utility-driven learning modification
term3 = learning_rate_term * learning_utility_term * \
        (gradient_of_loss + weight_regularization_term + decay_utility_term + external_stimulus_impact_term)

# Term 4: Utility-based pruning
term4 = pruning_probability_term * (-utility_based_pruning_term + pruning_randomness_term) # Simplified Activation(...)

# Term 5: Magnitude-based pruning
term5 = pruning_probability_term * weights_magnitude # Simplified Activation(|Weights|)

# Term 6: Global randomness
term6 = global_randomness_term * randomness_coefficient

# Term 7: Spatial diffusion
term7 = spatial_diffusion_term

# Term 8: Dynamic adaptive threshold
term8 = (base_threshold + fatigue_coefficient * recent_activity_integral - \
         cumulative_activity_coefficient * cumulative_activity_integral)

# Term 9: Memory trace
term9 = memory_decay_historical_influence_integral

# Term 10: Input relevance modulation
term10 = input_relevance_term * dropout_mask

# --- Calculate the final result for ∂w/∂t ---
dw_dt = term1 - term2 - term3 - term4 - term5 + term6 + term7 - term8 + term9 + term10

# --- Print the full equation with numerical values and the final result ---
print("Optimal Model: A")
print("Equation: ∂w(x, t) / ∂t = ")
print(f"  ({learning_rate_term:.2f} * ({mission_based_utility_term:.1f} + {gradient_of_loss:.1f}))")
print(f"- ({learning_rate_term:.2f} * ({gradient_of_loss:.1f} + {weight_regularization_term:.1f}))")
print(f"- ({learning_rate_term:.2f} * {learning_utility_term:.1f} * ({gradient_of_loss:.1f} + {weight_regularization_term:.1f} + {decay_utility_term:.1f} + {external_stimulus_impact_term:.1f}))")
print(f"- ({pruning_probability_term:.2f} * ({-utility_based_pruning_term:.1f} + {pruning_randomness_term:.1f}))")
print(f"- ({pruning_probability_term:.2f} * {weights_magnitude:.1f})")
print(f"+ ({global_randomness_term:.2f} * {randomness_coefficient:.1f})")
print(f"+ {spatial_diffusion_term:.1f}")
print(f"- ({base_threshold:.1f} + {fatigue_coefficient:.2f} * {recent_activity_integral:.1f} - {cumulative_activity_coefficient:.2f} * {cumulative_activity_integral:.1f})")
print(f"+ {memory_decay_historical_influence_integral:.2f}")
print(f"+ ({input_relevance_term:.1f} * {dropout_mask:.1f})")
print("\n--- Calculation ---")
print(f"∂w/∂t = ({term1:.4f}) - ({term2:.4f}) - ({term3:.4f}) - ({term4:.4f}) - ({term5:.4f}) + ({term6:.4f}) + ({term7:.4f}) - ({term8:.4f}) + ({term9:.4f}) + ({term10:.4f})")
print(f"∂w/∂t = {dw_dt:.4f}")
print("\nModel A is the optimal choice as it represents the most comprehensive and biologically plausible feature set for neuromorphic computing.")
