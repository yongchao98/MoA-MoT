import sys
import io

# Define a class to represent the primate's neurological state
class PrimateBrain:
    def __init__(self, lesion_side, lesion_pathway):
        self.lesion_side = lesion_side
        self.lesion_pathway = lesion_pathway
        self.conscious_deficit_quadrant = self._determine_deficit()

    def _determine_deficit(self):
        """Determines the visual field deficit based on the lesion."""
        if self.lesion_side == "right":
            field_side = "left"
        else:
            field_side = "right"

        if "meyer's loop" in self.lesion_pathway:
            field_vertical = "upper"
        elif "outside meyer's loop" in self.lesion_pathway:
            field_vertical = "lower"
        else:
            field_vertical = "unknown"
        
        return f"{field_vertical} {field_side}"

    def test_trial(self, stimulus_quadrant):
        """Simulates a trial and returns the primate's behavior."""
        is_in_deficit_field = (stimulus_quadrant == self.conscious_deficit_quadrant)
        
        # Blindsight means visuomotor pathways are intact despite lack of conscious vision
        if is_in_deficit_field:
            motor_action = f"Accurately reaches with {'left' if 'left' in stimulus_quadrant else 'right'} hand for the target."
            conscious_report = "Presses the button indicating 'no stimulus'."
            phenomenon = f"Blindsight for stimuli in the {self.conscious_deficit_quadrant} quadrant"
        else:
            motor_action = "Accurately reaches for the target."
            conscious_report = "Does not press the 'no stimulus' button."
            phenomenon = "Normal vision"
            
        return motor_action, conscious_report, phenomenon

# --- Main Simulation ---

# 1. Define the primate's condition based on the prompt
primate = PrimateBrain(lesion_side="right", lesion_pathway="white matter outside the meyer's loop")

# 2. Define the specific trial from the prompt
trial_stimulus_quadrant = "lower left"
trial_stimulus_details = f"A purple stimulus is given in the {trial_stimulus_quadrant} quadrant."

# 3. Run the simulation for this specific trial
motor_action, conscious_report, conclusion = primate.test_trial(trial_stimulus_quadrant)

# 4. Print the results clearly
print("--- Neuro-behavioral Simulation ---")
print(f"Condition: Ablation of {primate.lesion_side} {primate.lesion_pathway}.")
print(f"Predicted Conscious Deficit: Vision in the {primate.conscious_deficit_quadrant.title()} Quadrant is impaired.")
print("-" * 35)
print(f"Scenario: {trial_stimulus_details}")
print("\nObserved Behavior:")
print(f"  - Motor Action: {motor_action}")
print(f"  - 'Conscious' Report: {conscious_report}")
print("\nConclusion:")
print("The primate's ability to accurately act on a stimulus it does not consciously perceive is a demonstration of:")
print(f"'{conclusion}'")
