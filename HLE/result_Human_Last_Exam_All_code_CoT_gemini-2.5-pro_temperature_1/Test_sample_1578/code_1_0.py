import time

def simulate_mirror_test():
    """
    Simulates the test of using an opaque barrier to differentiate
    a true reflection from a demon mimic.
    """

    class TrueReflection:
        """Represents a real reflection, dependent on visibility."""
        def get_state(self, person_is_visible, person_action):
            if person_is_visible:
                return f"is visible and perfectly mimics your action, '{person_action}'."
            else:
                return "is obscured. The reflection is gone because you are hidden."

    class DemonMimic:
        """Represents a demon, an independent entity that mimics you."""
        def __init__(self):
            self.last_known_action = "standing still"

        def get_state(self, person_is_visible, person_action):
            if person_is_visible:
                self.last_known_action = person_action
                return f"is visible and mimics your action, '{person_action}'."
            else:
                # The demon can no longer see you and must guess or wait.
                return (f"remains visible, frozen in your last action, '{self.last_known_action}'. "
                        "This is not a reflection!")

    # Setup
    reflection = TrueReflection()
    demon = DemonMimic()
    person_action = "waving your hand"
    person_is_visible = True

    print("You stand before the mirror.")
    print(f"You are currently: {person_action}\n")

    print("--- Testing Scenario 1: It is a true reflection ---")
    print(f"The figure in the mirror {reflection.get_state(person_is_visible, person_action)}")
    print("\n--- Testing Scenario 2: It is a demon ---")
    print(f"The figure in the mirror {demon.get_state(person_is_visible, person_action)}")

    print("\n" + "="*50 + "\n")
    print("You now raise an opaque barrier, hiding yourself completely.")
    print("Behind the barrier, you are now 'patting your head'.")
    
    # Update state after raising the barrier
    person_is_visible = False
    person_action = "patting your head" # The entity in the mirror cannot know this.

    print("\n--- Result for Scenario 1 (True Reflection) ---")
    print(f"The figure in the mirror {reflection.get_state(person_is_visible, person_action)}")

    print("\n--- Result for Scenario 2 (Demon) ---")
    print(f"The figure in the mirror {demon.get_state(person_is_visible, person_action)}")

    print("\n\nConclusion: By blocking the line of sight, the demon's inability to "
          "replicate an unseen action reveals its true nature.")


simulate_mirror_test()