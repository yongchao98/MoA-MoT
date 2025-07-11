import sys

class PrimateVisualSystem:
    """
    A class to simulate the visual system of a primate and the effects of a lesion.
    """
    def __init__(self):
        """
        Initializes the visual field quadrants.
        'aware' means conscious perception is intact.
        """
        self.visual_field = {
            'upper_right': 'aware',
            'upper_left': 'aware',
            'lower_right': 'aware',
            'lower_left': 'aware'
        }
        print("Primate visual system initialized. All quadrants are perceived consciously.")
        print("-" * 20)

    def perform_lesion(self):
        """
        Simulates the described radiofrequency thermal ablation.
        This lesion affects the pathway for the lower-left visual field.
        """
        lesion_area = "right optic radiation (non-Meyer's loop fibers)"
        affected_quadrant = "lower_left"
        
        print(f"Performing ablation on: {lesion_area}")
        print(f"This pathway processes the contralateral inferior visual field.")
        print(f"Result: Conscious visual pathway for the '{affected_quadrant}' quadrant is disrupted.")
        
        # Update the state of the affected quadrant to 'blindsight'
        self.visual_field[affected_quadrant] = 'blindsight'
        print(f"The state of the '{affected_quadrant}' quadrant is now '{self.visual_field[affected_quadrant]}'.")
        print("-" * 20)

    def test_stimulus(self, location, color):
        """
        Tests the primate's response to a visual stimulus in a specific location.
        """
        print(f"A {color} stimulus is presented in the {location.replace('_', ' ')} quadrant.")
        
        quadrant_state = self.visual_field.get(location)

        if quadrant_state == 'blindsight':
            print("\nDemonstrated Behavior:")
            print("1. Action: The primate accurately reaches with its left hand for the target.")
            print("2. Report: The primate then presses the button indicating 'no trial' or 'no stimulus seen'.")
            print("\nConclusion:")
            print("This demonstrates blindsight for stimuli in the lower left quadrant in a non-verbal primate.")
            
        elif quadrant_state == 'aware':
            print("\nDemonstrated Behavior:")
            print("1. Action: The primate accurately reaches for the target.")
            print("2. Report: The primate indicates it consciously saw the stimulus.")
            print("\nConclusion:")
            print(f"This demonstrates normal vision in the {location.replace('_', ' ')} quadrant.")
            
        else:
            print("Invalid stimulus location.")

# --- Main execution ---
if __name__ == "__main__":
    # Create an instance of the primate's visual system
    primate_brain = PrimateVisualSystem()

    # Perform the lesion as described in the problem
    primate_brain.perform_lesion()

    # Test the specific scenario from the question
    stimulus_location = 'lower_left'
    stimulus_color = 'purple'
    primate_brain.test_stimulus(stimulus_location, stimulus_color)
