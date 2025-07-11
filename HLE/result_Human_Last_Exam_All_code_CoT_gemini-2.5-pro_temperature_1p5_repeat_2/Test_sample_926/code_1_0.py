import sys
import io

# Buffer to capture the output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def calculate_friction(normal_load, coherence_factor, base_coefficient=0.5):
    """
    Calculates the frictional force based on a simplified model of superlubricity.

    Args:
        normal_load (float): The force pressing the surfaces together (in nN).
        coherence_factor (float): A value representing atomic lattice alignment.
                                  - 1.0 for high coherence (aligned, high friction).
                                  - A small value (e.g., 0.01) for low coherence
                                    (misaligned, superlubricity).
        base_coefficient (float): A base friction coefficient.

    Returns:
        float: The calculated frictional force (in nN).
    """
    friction = base_coefficient * normal_load * coherence_factor
    return friction

def demonstrate_superlubricity_principles():
    """
    Demonstrates how normal load and coherence affect friction in superlubric systems.
    """
    print("This model demonstrates the principles of superlubricity from option B.")
    print("Frictional Force = Base Coefficient * Normal Load * Coherence Factor\n")

    # --- Scenario Parameters ---
    load_low = 10.0  # nanoNewtons
    load_high = 50.0 # nanoNewtons
    coeff = 0.5
    coherence_high = 1.0  # Represents aligned atomic lattices (high friction)
    coherence_low = 0.01 # Represents misaligned lattices (superlubricity)

    # --- Calculations & Explanations ---

    # Scenario 1: Low load, but high coherence (no superlubricity)
    friction1 = calculate_friction(load_low, coherence_high, coeff)
    print(f"1. High Coherence Case (Aligned Lattices):")
    print(f"   Friction = {coeff} * {load_low} nN * {coherence_high} = {friction1:.2f} nN")
    print("   Result: Even with a low load, high atomic coherence leads to significant friction.\n")


    # Scenario 2: Low load and low coherence (superlubricity)
    friction2 = calculate_friction(load_low, coherence_low, coeff)
    print(f"2. Low Coherence Case (Superlubricity):")
    print(f"   Friction = {coeff} * {load_low} nN * {coherence_low} = {friction2:.2f} nN")
    print("   Result: Low coherence (misalignment) drastically reduces friction, enabling superlubricity.\n")


    # Scenario 3: Higher load but still low coherence (superlubric state under more pressure)
    friction3 = calculate_friction(load_high, coherence_low, coeff)
    print(f"3. Superlubricity with Higher Load:")
    print(f"   Friction = {coeff} * {load_high} nN * {coherence_low} = {friction3:.2f} nN")
    print("   Result: Increasing the normal load increases the friction, even in the superlubric state.\n")

    print("Conclusion: Frictional force increases with BOTH normal load and atomic coherence.")


# Run the demonstration
demonstrate_superlubricity_principles()

# Restore stdout and print the captured output
sys.stdout = old_stdout
print(captured_output.getvalue())
<<<B>>>