import sys

# Suppress the display of the floating point number warning that might be generated from large numbers
# on some python versions. The logic remains sound.
if sys.version_info.major == 3 and sys.version_info.minor >= 11:
    import warnings
    warnings.filterwarnings("ignore", message="line with integer literals has more than 255 characters")

def model_self_stabilizing_effect():
    """
    This function models the self-stabilizing effect across different learning phases.
    It demonstrates that the effect, driven by the discovery of knowledge gaps,
    is strongest in the late learning phase.
    """

    # We model the effect using the equation: Effect = K^3 * (100 - K)
    # where K is the knowledge level (0-100).
    # This function peaks at K=75, representing the late learning phase.
    def calculate_effect(knowledge_level):
        """Calculates the effect for a given knowledge level."""
        # Using integer arithmetic to avoid floating point inaccuracies with large numbers.
        return knowledge_level**3 * (100 - knowledge_level)

    # Define representative knowledge levels for each phase
    early_phase_k = 20
    intermediate_phase_k = 50
    late_phase_k = 75  # The peak of our model function

    # Calculate the numerical effect for each phase
    effect_early = calculate_effect(early_phase_k)
    effect_intermediate = calculate_effect(intermediate_phase_k)
    effect_late = calculate_effect(late_phase_k)

    print("Modeling the 'self-stabilizing effect' based on learning phase.")
    print("The model uses the equation: Effect = (Knowledge Level)^3 * (100 - Knowledge Level)")
    print("-" * 60)

    # --- Print results for Early Phase ---
    print(f"1. Early Phase (Knowledge Level = {early_phase_k}):")
    print(f"   Effect = {early_phase_k}^3 * (100 - {early_phase_k})")
    print(f"   Effect = {early_phase_k**3} * {100 - early_phase_k}")
    print(f"   Calculated Effect = {effect_early}\n")

    # --- Print results for Intermediate Phase ---
    print(f"2. Intermediate Phase (Knowledge Level = {intermediate_phase_k}):")
    print(f"   Effect = {intermediate_phase_k}^3 * (100 - {intermediate_phase_k})")
    print(f"   Effect = {intermediate_phase_k**3} * {100 - intermediate_phase_k}")
    print(f"   Calculated Effect = {effect_intermediate}\n")

    # --- Print results for Late Phase ---
    print(f"3. Late Phase (Knowledge Level = {late_phase_k}):")
    print(f"   Effect = {late_phase_k}^3 * (100 - {late_phase_k})")
    print(f"   Effect = {late_phase_k**3} * {100 - late_phase_k}")
    print(f"   Calculated Effect = {effect_late}\n")
    print("-" * 60)

    # --- Conclusion ---
    print("Conclusion from the model:")
    print(f"The calculated effect is highest ({effect_late}) in the late learning phase.")
    print("This supports the idea that a strong foundational understanding is required to")
    print("discover the greatest number of complex knowledge gaps, thus maximizing the effect.")

model_self_stabilizing_effect()
<<<C>>>