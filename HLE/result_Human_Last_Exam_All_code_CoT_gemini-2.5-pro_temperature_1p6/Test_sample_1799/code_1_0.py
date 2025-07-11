def demonstrate_postulates():
    """
    This function demonstrates why the second postulate of special relativity
    is not considered superfluous.
    """
    # Speed of light in a vacuum (m/s)
    c = 299792458

    # Let's imagine a spaceship traveling at a significant fraction of the speed of light
    # e.g., 20% of the speed of light
    v_spaceship = 0.20 * c

    print("Imagine a spaceship firing a laser beam in its direction of travel.")
    print(f"The spaceship is moving at v = {int(v_spaceship)} m/s (0.20c).")
    print("An observer on a stationary platform measures the speed of the light from the laser.")
    print("-" * 70)

    # --- Scenario based on Postulate 1 alone (assuming classical physics) ---
    print("Scenario 1: Using only the Principle of Relativity (Postulate 1)")
    print("This allows for Galilean relativity, where velocities add up.")
    print("The predicted speed of light for the stationary observer would be c + v.")
    
    # Calculate the Galilean prediction
    c_observed_galilean = c + v_spaceship
    
    print("\nEquation (Galilean): Observed Speed = c + v")
    print(f"Predicted Speed = {c} + {int(v_spaceship)}")
    print(f"Predicted Speed = {int(c_observed_galilean)} m/s")
    print("-" * 70)

    # --- Scenario based on Postulates 1 and 2 (Special Relativity) ---
    print("Scenario 2: Using both Postulate 1 and Postulate 2")
    print("The 2nd postulate states the speed of light is constant for all observers.")
    print("The predicted speed of light for the stationary observer would be c.")
    
    # The relativistic prediction is simply c
    c_observed_relativistic = c

    print("\nEquation (Relativistic): Observed Speed = c")
    print(f"Predicted Speed = {c}")
    print(f"Predicted Speed = {c_observed_relativistic} m/s")
    print("-" * 70)
    
    print("\nConclusion:")
    print("The two postulates lead to different, testable predictions. Experiments have")
    print("consistently confirmed Scenario 2. Therefore, the second postulate is a")
    print("necessary, independent assertion and is not superfluous.")

if __name__ == '__main__':
    demonstrate_postulates()