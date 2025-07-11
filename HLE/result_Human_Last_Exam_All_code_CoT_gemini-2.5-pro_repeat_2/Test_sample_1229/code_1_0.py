# --- Assigning example numerical values to each term ---

# Learning and Gradient Terms
learning_rate_term = 0.01
mission_based_utility_term = 5.0
gradient_of_loss_with_respect_to_weights = 2.0
weight_regularization_term = 0.1

# Utility-Modulated Learning Terms
learning_utility_term = 0.5
decay_utility_term = 0.2
external_stimulus_impact_term = 0.3

# Pruning Terms
pruning_probability_term = 0.05
# For simplicity, we assume ActivationFunction(x) = x
# Term 1: Activation Function(− Utility-Based Pruning Term + Randomness Term)
utility_based_pruning_term = 1.0
randomness_term_pruning = 0.4
# Term 2: Activation Function (|Weights|)
weights_magnitude = 0.8

# Stochastic and Spatial Terms
global_randomness_term = 0.001
randomness_coefficient = 1.5
spatial_diffusion_term = 0.02

# Adaptive Threshold Terms (integrals are represented as pre-calculated values)
base_threshold = 0.05
fatigue_coefficient = 0.2
integral_recent_activity = 0.6
cumulative_activity_coefficient = 0.01
integral_cumulative_activity = 10.0

# Memory and Input Terms
integral_memory_decay = 0.03
input_relevance_x_dropout = 0.015


# --- Calculate each major component of Equation A ---

term1 = learning_rate_term * (mission_based_utility_term + gradient_of_loss_with_respect_to_weights)
term2 = -learning_rate_term * (gradient_of_loss_with_respect_to_weights + weight_regularization_term)
term3 = -learning_rate_term * learning_utility_term * (gradient_of_loss_with_respect_to_weights + weight_regularization_term + decay_utility_term + external_stimulus_impact_term)
term4 = -pruning_probability_term * (-utility_based_pruning_term + randomness_term_pruning)
term5 = -pruning_probability_term * weights_magnitude
term6 = global_randomness_term * randomness_coefficient
term7 = spatial_diffusion_term
term8 = -(base_threshold + fatigue_coefficient * integral_recent_activity - cumulative_activity_coefficient * integral_cumulative_activity)
term9 = integral_memory_decay
term10 = input_relevance_x_dropout

# Calculate the final result for Differential Updates (∂w(x, t) / ∂t)
differential_update = term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8 + term9 + term10

# --- Print the final equation with all numbers substituted ---

print("Model A: Differential Updates ( ∂w(x, t) / ∂t ) = ")
print(f"  {learning_rate_term} * ({mission_based_utility_term} + {gradient_of_loss_with_respect_to_weights})")
print(f"  - {learning_rate_term} * ({gradient_of_loss_with_respect_to_weights} + {weight_regularization_term})")
print(f"  - {learning_rate_term} * {learning_utility_term} * ({gradient_of_loss_with_respect_to_weights} + {weight_regularization_term} + {decay_utility_term} + {external_stimulus_impact_term})")
print(f"  - {pruning_probability_term} * ActivationFunction(-{utility_based_pruning_term} + {randomness_term_pruning})")
print(f"  - {pruning_probability_term} * ActivationFunction(|{weights_magnitude}|)")
print(f"  + {global_randomness_term} * {randomness_coefficient}")
print(f"  + {spatial_diffusion_term}")
print(f"  - ({base_threshold} + {fatigue_coefficient} * {integral_recent_activity} - {cumulative_activity_coefficient} * {integral_cumulative_activity})")
print(f"  + {integral_memory_decay}")
print(f"  + {input_relevance_x_dropout}")
print("\nCalculating the result with the given numbers...")
print(f"Differential Update Value = {differential_update:.6f}")
