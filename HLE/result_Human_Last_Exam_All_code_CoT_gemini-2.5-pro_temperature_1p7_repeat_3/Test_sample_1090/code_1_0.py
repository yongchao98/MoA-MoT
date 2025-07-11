def check_orthostatic_hypotension():
    """
    Analyzes patient vital signs to check for orthostatic hypotension.
    """
    # Patient's vital signs
    resting_systolic_bp = 135
    resting_diastolic_bp = 85
    standing_systolic_bp = 150
    standing_diastolic_bp = 90

    # Diagnostic criteria for orthostatic hypotension (a drop in BP)
    systolic_drop_threshold = 20
    diastolic_drop_threshold = 10

    # Calculate the change in blood pressure.
    # A positive value indicates a drop in pressure.
    systolic_change = resting_systolic_bp - standing_systolic_bp
    diastolic_change = resting_diastolic_bp - standing_diastolic_bp

    print("Step 1: Analyzing for Orthostatic Hypotension")
    print("="*45)
    print(f"Patient's resting blood pressure: {resting_systolic_bp}/{resting_diastolic_bp} mmHg")
    print(f"Patient's standing blood pressure: {standing_systolic_bp}/{standing_diastolic_bp} mmHg")
    print("-" * 45)

    # Print the equation for systolic change
    print(f"Equation for systolic pressure change (drop):")
    print(f"{resting_systolic_bp} (Resting) - {standing_systolic_bp} (Standing) = {systolic_change} mmHg")

    # Print the equation for diastolic change
    print(f"\nEquation for diastolic pressure change (drop):")
    print(f"{resting_diastolic_bp} (Resting) - {standing_diastolic_bp} (Standing) = {diastolic_change} mmHg")
    print("-" * 45)

    # Check if the criteria for orthostatic hypotension are met
    is_systolic_orthostatic = systolic_change >= systolic_drop_threshold
    is_diastolic_orthostatic = diastolic_change >= diastolic_drop_threshold

    print("Step 2: Conclusion from Vitals Analysis")
    print("="*45)
    if is_systolic_orthostatic or is_diastolic_orthostatic:
        print("Result: The patient HAS orthostatic hypotension.")
    else:
        print("Result: The patient does NOT have orthostatic hypotension.")
        print(f"The systolic pressure did not drop; it changed by {systolic_change} mmHg.")
        print(f"The diastolic pressure did not drop; it changed by {diastolic_change} mmHg.")
        print("\nClinical Reasoning: Since orthostatic hypotension is ruled out, the inability to ambulate is likely due to another cause. The significant resistance to knee extension points to post-stroke spasticity as the primary limiting factor. The most appropriate next step is to address this directly to enable effective physical therapy.")

# Run the analysis
check_orthostatic_hypotension()