def find_diagnostic_maneuver():
    """
    Analyzes the clinical vignette to determine the most likely diagnosis
    and the physical exam maneuver that would confirm it.
    """

    # --- Step 1: Analyze Patient Presentation ---
    print("Step 1: Analyzing the patient's key clinical features.")
    pain_location = "L4-S1"
    pain_exacerbator = "lying supine"
    imaging_result = "unremarkable X-ray"
    
    print(f"The patient has pain in the {pain_location} distribution, which suggests sciatic nerve involvement.")
    print(f"The {imaging_result} result makes a primary bony pathology less likely, pointing towards a soft-tissue cause.")
    print("-" * 30)

    # --- Step 2: Formulate the Most Likely Diagnosis ---
    print("Step 2: Formulating the diagnostic hypothesis.")
    diagnosis = "Piriformis Syndrome"
    reasoning = "This condition occurs when the piriformis muscle irritates the sciatic nerve, mimicking radicular pain without spinal pathology."
    print(f"The most likely diagnosis fitting the symptoms and normal X-ray is {diagnosis}.")
    print(f"Reasoning: {reasoning}")
    print("-" * 30)

    # --- Step 3: Identify the Confirmatory Test ---
    print("Step 3: Determining the confirmatory physical exam test.")
    muscle_to_test = "piriformis muscle"
    primary_action = "External Rotation of the hip"
    
    print(f"To confirm {diagnosis}, the test must stress the {muscle_to_test}.")
    print(f"The primary action of the piriformis muscle is: {primary_action}.")
    print("Therefore, applying resistance while the patient attempts this action would reproduce the pain if the piriformis muscle is the cause.")
    print("-" * 30)
    
    # --- Step 4: Final Conclusion ---
    answer_option = "D"
    answer_description = "External Rotation"
    print("Step 4: Conclusion.")
    print(f"The action that will confirm the diagnosis is '{answer_description}'.")
    print(f"This corresponds to answer choice {answer_option}.")

# Run the analysis
find_diagnostic_maneuver()