def diagnose_murmur():
    """
    This script analyzes the clinical findings from the case to determine the most likely cause of the heart murmur.
    It will print its reasoning step-by-step.
    """
    # Step 1: Define the key findings from the clinical vignette.
    print("Step 1: Listing the key clinical findings.")
    murmur_location = "Left upper sternal border (pulmonic area)"
    murmur_type = "Systolic ejection murmur"
    dynamic_maneuver = "Increases in intensity with inspiration"
    ecg_finding = "Right Ventricular Hypertrophy (RVH) with Left Axis Deviation (LAD)"
    
    # Step 2: Create a 'diagnostic equation' by interpreting the findings.
    print("\nStep 2: Building the diagnostic equation from the findings.")
    print("Finding 1: The murmur location and type point to the pulmonic valve.")
    print(f"Equation Part 1: ( {murmur_location} ) + ( {murmur_type} ) => Murmur from Pulmonic Valve")
    
    print("\nFinding 2: The change with breathing confirms a right-sided heart problem.")
    print(f"Equation Part 2: ( Murmur from Pulmonic Valve ) + ( {dynamic_maneuver} ) => Right-Sided Heart Pathology")

    print("\nFinding 3: Evaluating the cause of the pulmonic murmur.")
    print("This murmur is due to high blood flow across the pulmonic valve (relative stenosis).")
    print("Among the choices, an Atrial Septal Defect (ASD) causes a large left-to-right shunt, creating this exact high-flow state.")
    print("This also explains the Right Ventricular Hypertrophy (RVH).")
    
    # Step 3: Conclude the diagnosis by assembling the equation.
    print("\nStep 3: Final Diagnostic Equation.")
    print(f"( {murmur_type} at {murmur_location} ) + ( {dynamic_maneuver} ) + ( {ecg_finding} ) => Atrial Septal Defect")
    print("\nExplanation: The combination of findings, especially a flow murmur across the pulmonic valve and the specific ECG pattern (LAD + RVH), is highly suggestive of an ostium primum Atrial Septal Defect.")

diagnose_murmur()
<<<D>>>