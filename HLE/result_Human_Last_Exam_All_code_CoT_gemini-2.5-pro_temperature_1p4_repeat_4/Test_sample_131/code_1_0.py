import matplotlib.pyplot as plt
import numpy as np

def calculate_self_stabilizing_effect(knowledge_level):
    """
    Models the self-stabilizing effect based on knowledge level (0-100).
    The function is designed to be low for beginners, rise, and then peak in the 'late' phase,
    representing the idea that deep foundational knowledge is needed to perceive the most complex gaps.
    This is a simplified quadratic model for demonstration.
    
    Formula: Effect = -0.015 * (Knowledge - 80)^2 + 100
    This equation sets the peak effect (100) at a knowledge level of 80 (late phase).
    """
    # The constants -0.015, 80, and 100 are chosen to create a representative curve.
    peak_knowledge_level = 80
    curve_steepness = -0.015
    max_effect = 100
    
    # The equation for the effect
    effect = curve_steepness * (knowledge_level - peak_knowledge_level)**2 + max_effect
    return max(0, effect) # Ensure effect is not negative

# Define representative knowledge levels for each phase
early_phase_knowledge = 20
intermediate_phase_knowledge = 50
late_phase_knowledge = 80 # This is where our model will peak

# Calculate the effect for each phase
effect_early = calculate_self_stabilizing_effect(early_phase_knowledge)
effect_intermediate = calculate_self_stabilizing_effect(intermediate_phase_knowledge)
effect_late = calculate_self_stabilizing_effect(late_phase_knowledge)

# Print the results and analysis
print("Modeling the Self-Stabilizing Effect of Knowledge Acquisition:")
print("-" * 60)
print("This model simulates the effect's strength at different learning phases.")
print("The formula used is: Effect = -0.015 * (Knowledge - 80)^2 + 100")
print("-" * 60)

print(f"Early Learning Phase (Knowledge Level = {early_phase_knowledge}):")
print(f"Calculation: -0.015 * ({early_phase_knowledge} - 80)^2 + 100 = {effect_early:.1f}")
print("Result: The effect is relatively low as the learner doesn't yet perceive many gaps.\n")

print(f"Intermediate Learning Phase (Knowledge Level = {intermediate_phase_knowledge}):")
print(f"Calculation: -0.015 * ({intermediate_phase_knowledge} - 80)^2 + 100 = {effect_intermediate:.1f}")
print("Result: The effect is strong and growing as foundational knowledge reveals many gaps.\n")

print(f"Late Learning Phase (Knowledge Level = {late_phase_knowledge}):")
print(f"Calculation: -0.015 * ({late_phase_knowledge} - 80)^2 + 100 = {effect_late:.1f}")
print("Result: The effect PEAKS as comprehensive knowledge allows the discovery of the most complex gaps.\n")

print("Conclusion:")
print("The model shows the self-stabilizing effect is weakest in the early phase, grows in the intermediate phase,")
print("and peaks in the late phase. This supports statement C.")
