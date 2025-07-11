def model_superlubricity_friction(normal_load, alignment_factor):
    """
    Models the frictional force in a system exhibiting superlubricity.

    Args:
        normal_load (float): The load pressing the surfaces together (in nanoNewtons, nN).
        alignment_factor (float): A value from 0.0 (perfectly misaligned, incommensurate)
                                  to 1.0 (perfectly aligned, commensurate).
    """
    # A base friction coefficient for the ideal superlubric state.
    base_coefficient = 0.001

    # The model shows that friction is proportional to the normal load.
    # The alignment factor introduces a multiplier effect. When alignment is 0,
    # the multiplier is 1 (ideal superlubricity). When alignment is 1,
    # the multiplier is large, representing high friction in a commensurate state.
    friction_multiplier = 1 + alignment_factor * 100
    frictional_force = base_coefficient * normal_load * friction_multiplier

    print(f"--- Calculating Friction ---")
    print(f"Inputs: Normal Load = {normal_load} nN, Alignment Factor = {alignment_factor}")
    print("Equation: Frictional Force = Base Coefficient * Normal Load * (1 + Alignment Factor * 100)")
    # Final equation with numbers plugged in
    print(f"Calculation: {frictional_force:.4f} nN = {base_coefficient} * {normal_load} * (1 + {alignment_factor} * 100)")
    print("-" * 28 + "\n")


# Case 1: Superlubric State (low alignment, low load)
# Friction is extremely low.
model_superlubricity_friction(normal_load=5.0, alignment_factor=0.0)

# Case 2: Increasing the Normal Load
# With alignment still low, friction increases proportionally to the load.
model_superlubricity_friction(normal_load=20.0, alignment_factor=0.0)

# Case 3: Increasing the Alignment
# Even with a low load, increasing alignment dramatically increases friction.
model_superlubricity_friction(normal_load=5.0, alignment_factor=0.8)

# Case 4: High Friction State (high alignment, high load)
# Both factors contribute to a much higher frictional force.
model_superlubricity_friction(normal_load=20.0, alignment_factor=1.0)