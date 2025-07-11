import sys
# Redirect print to a string buffer to check if it's too long
from io import StringIO
original_stdout = sys.stdout
sys.stdout = captured_output = StringIO()

def solve():
    """
    This script identifies the dmPFC sub-region most purely activated by emotion processing.
    It does so by modeling the known functional specializations of each region and
    calculating an "emotion purity" score.
    """
    
    # Step 1: Model the functional profiles based on neuroscientific literature.
    # Scores are illustrative, representing relative engagement in different domains.
    # A high score in 'Emotion' and low scores elsewhere signifies high purity.
    dmPFC_regions = {
        'rostroventral': {
            'Description': 'Processes value, affect, and self-related emotion.',
            'Emotion': 75,
            'Cognition': 20,
            'Motor': 5
        },
        'rostrodorsal': {
            'Description': 'Involved in metacognition and cognitive appraisal of emotion.',
            'Emotion': 60,
            'Cognition': 90,
            'Motor': 10
        },
        'caudal-left': {
            'Description': 'Related to higher-order motor control and action selection.',
            'Emotion': 25,
            'Cognition': 45,
            'Motor': 75
        },
        'caudal-right': {
            'Description': 'Related to higher-order motor control and action selection.',
            'Emotion': 20,
            'Cognition': 40,
            'Motor': 80
        }
    }

    # Step 2 & 3: Define a dictionary to hold results and calculate purity for each region.
    results = {}
    print("Analysis of Emotion Processing Purity in dmPFC Sub-regions\n" + "="*60)
    
    for region, profile in dmPFC_regions.items():
        emotion_act = profile['Emotion']
        cognition_act = profile['Cognition']
        motor_act = profile['Motor']
        
        # Calculate the total functional engagement
        total_act = emotion_act + cognition_act + motor_act
        
        # Calculate emotion purity
        purity = emotion_act / total_act if total_act > 0 else 0
        
        results[region] = purity
        
        # Print the step-by-step calculation for transparency
        print(f"Calculating for region: '{region}'")
        print(f"Profile: {profile['Description']}")
        print(f"Purity Formula: Emotion / (Emotion + Cognition + Motor)")
        print(f"Calculation: {emotion_act} / ({emotion_act} + {cognition_act} + {motor_act}) = {purity:.4f}\n")

    # Step 4: Identify the region with the highest purity score.
    most_pure_region = max(results, key=results.get)
    max_purity_score = results[most_pure_region]
    final_calculation_profile = dmPFC_regions[most_pure_region]

    print("="*60 + "\nConclusion:\n")
    print(f"The region most purely activated by emotion processing is '{most_pure_region}'.")
    print(f"It has the highest emotion purity score of {max_purity_score:.4f}, based on the following functional profile calculation:")
    
    # Fulfilling the requirement to show each number in the final equation again.
    e_val = final_calculation_profile['Emotion']
    c_val = final_calculation_profile['Cognition']
    m_val = final_calculation_profile['Motor']
    print(f"Final Equation: {e_val} / ({e_val} + {c_val} + {m_val})")

solve()

# Restore original stdout and print the captured output
sys.stdout = original_stdout
output_str = captured_output.getvalue()
print(output_str)
