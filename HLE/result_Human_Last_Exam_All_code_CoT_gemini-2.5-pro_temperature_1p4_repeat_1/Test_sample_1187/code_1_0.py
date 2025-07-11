import sys

# Define known dermatome landmarks as a dictionary for reference.
# A dermatome is an area of skin that is mainly supplied by a single spinal nerve.
dermatome_landmarks = {
    'nipple_line': 'T4',
    'xiphoid_process': 'T6',
    'umbilicus': 'T10',
    'groin': 'T12'
}

# Step 1: Analyze the patient's neurological symptoms from the vignette.
patient_symptoms = {
    'right_leg_motor': 'severe weakness (corticospinal tract)',
    'right_leg_sensation': 'loss of proprioception and vibration (dorsal columns)',
    'left_side_sensation': 'loss of pain and temperature (spinothalamic tract)',
    'sensory_level': 'umbilicus'
}

# Step 2: Identify the key localizing sign.
# The most reliable indicator for the vertical level of a spinal cord injury
# is the sensory level, which is the highest point of sensory loss.
key_localizing_sign = patient_symptoms['sensory_level']

print(f"Analyzing the clinical presentation:")
print(f"- The pattern of ipsilateral motor/vibration loss and contralateral pain/temperature loss points to a spinal cord hemisection (Brown-Séquard Syndrome).")
print(f"- The key to determining the location is the sensory level.")
print(f"- The sensory loss begins at the '{key_localizing_sign}'.")
print("-" * 20)


# Step 3: Map the sensory level to the corresponding spinal cord segment.
if key_localizing_sign in dermatome_landmarks:
    lesion_level = dermatome_landmarks[key_localizing_sign]
    lesion_number_str = ''.join(filter(str.isdigit, lesion_level))
    lesion_number = int(lesion_number_str)
    
    print(f"The sensory level at the umbilicus corresponds to the dermatome T10.")
    print(f"This indicates the spinal cord injury is at the T10 level.")
    print("-" * 20)
    
    # Step 4: Present the final conclusion in the requested format.
    print(f"Final Equation Output:")
    print(f"Sensory Level ({key_localizing_sign}) = Lesion Level (T{lesion_number})")
    
else:
    print("Could not determine the lesion level based on the provided sensory level.", file=sys.stderr)
