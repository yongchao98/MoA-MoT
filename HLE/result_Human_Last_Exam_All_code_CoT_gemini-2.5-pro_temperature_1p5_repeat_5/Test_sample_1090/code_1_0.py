# This script simulates the clinical reasoning to determine the next step in management.

# --- Patient Data ---
# Key information from the clinical vignette.
history_of_stroke = True
recent_prolonged_immobility = True
is_unable_to_ambulate_despite_pt = True

# The most critical new physical exam finding.
full_knee_extension_degrees = 180
patient_knee_extension_limit_degrees = 165

# --- Analysis ---
# The primary reason for the patient's inability to ambulate is a new physical limitation.
# We will quantify this limitation to guide our decision.
print("Step 1: Quantify the new physical deficit.")
knee_extension_deficit = full_knee_extension_degrees - patient_knee_extension_limit_degrees

# Per the user request, showing the numbers in the final equation.
print(f"A fully extended knee is {full_knee_extension_degrees} degrees.")
print(f"The patient's knee extension is limited to {patient_knee_extension_limit_degrees} degrees.")
print(f"Calculating the extension deficit: {full_knee_extension_degrees} - {patient_knee_extension_limit_degrees} = {knee_extension_deficit} degrees.")

# --- Conclusion and Recommendation ---
print("\nStep 2: Determine the most appropriate action based on the findings.")
# The combination of stroke history, immobility, and a new contracture points to a specific problem.
if knee_extension_deficit > 0 and history_of_stroke and recent_prolonged_immobility:
    print("\nClinical Conclusion:")
    print("The patient's inability to walk is most directly caused by a developing knee contracture due to worsened post-stroke spasticity after a period of bed rest.")
    
    # This specific problem requires specialized management.
    next_step = "Physical medicine and rehabilitation consultation"
    print("\nSingle Most Appropriate Next Step in Management:")
    print(next_step)
else:
    # This case is clear, but a fallback is included.
    next_step = "Further non-specific investigation is required."
    print(next_step)
