import sys

def find_purest_emotion_region():
    """
    This script simulates a neuroimaging meta-analysis to find the sub-region of the
    dorsomedial prefrontal cortex (dmPFC) most purely activated by emotion.

    "Purity" is defined as the activation for emotion relative to activation for
    other major functions (cognition and motor control).

    The data is hypothetical but reflects established functional specializations:
    - rostroventral: Core emotion/valuation.
    - rostrodorsal: Social cognition, theory of mind (has emotional component but is less pure).
    - caudal: Motor planning/control (pre-SMA/SMA).
    """

    # Hypothetical activation data for dmPFC sub-regions.
    # Note: 'rostrodorsal' has the highest absolute emotion value (80),
    # but 'rostroventral' is more selective for emotion.
    activations = {
        "caudal-right": {"Emotion": 20, "Cognition": 30, "Motor": 80},
        "rostroventral": {"Emotion": 70, "Cognition": 20, "Motor": 5},
        "rostrodorsal": {"Emotion": 80, "Cognition": 85, "Motor": 10},
        "caudal-left": {"Emotion": 22, "Cognition": 28, "Motor": 85},
    }

    max_purity = -1
    purest_region = None
    final_equation_vars = {}

    print("Calculating emotion processing purity for each dmPFC sub-region...")
    print("Formula: Purity = Emotion / (Cognition + Motor)\n")

    for region, tasks in activations.items():
        emotion_act = tasks["Emotion"]
        cognition_act = tasks["Cognition"]
        motor_act = tasks["Motor"]

        # To avoid division by zero, although not expected with this data
        other_activation_sum = cognition_act + motor_act
        if other_activation_sum == 0:
            # In the unlikely event of zero other activations, purity is effectively infinite
            purity = float('inf')
        else:
            purity = emotion_act / other_activation_sum

        print(f"- {region}: Purity = {purity:.2f}")

        if purity > max_purity:
            max_purity = purity
            purest_region = region
            final_equation_vars = {
                "Emotion": emotion_act,
                "Cognition": cognition_act,
                "Motor": motor_act
            }

    print("\n---")
    print(f"The sub-region most purely activated by emotion processing is: {purest_region}")
    
    # Outputting the final equation as requested
    emotion_val = final_equation_vars["Emotion"]
    cognition_val = final_equation_vars["Cognition"]
    motor_val = final_equation_vars["Motor"]
    
    print(f"This was determined by finding the maximum purity score, calculated for the winning region as follows:")
    print(f"Purity = {emotion_val} / ({cognition_val} + {motor_val}) = {max_purity:.2f}")

# Execute the function
find_purest_emotion_region()
<<<B>>>