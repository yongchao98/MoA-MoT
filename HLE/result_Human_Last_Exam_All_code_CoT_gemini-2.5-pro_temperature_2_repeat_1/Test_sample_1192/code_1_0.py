import textwrap

def analyze_post_op_care():
    """
    Analyzes a clinical scenario to determine the best course of action
    to prevent post-operative complications after heart valve surgery.
    """
    
    # The clinical scenario details
    scenario = "A 56-year-old man is doing extremely well after undergoing heart valve surgery. His vitals are stable, and he feels fine. He wants to go home."
    
    # The question to be answered
    question = "What is the next course of action to prevent adverse post-operative complications?"

    # Define the answer choices
    choices = {
        'A': "Do not prescribe any medication since patient is asymptomatic and doing well.",
        'B': "Prescribe an analgesic for breakthrough pain.",
        'C': "Schedule physical therapy appointment.",
        'D': "Encourage regular exercise to increase circulation.",
        'E': "Return to the clinic in one month for routine post-operative follow-up.",
        'F': "Since only symptomatic treatment is recommended, no action is needed at this time.",
        'G': "Keep patient at the hospital for one more day.",
        'H': "Discharge patient with heart-healthy dietary instructions.",
        'I': "None of the answer choices.",
        'J': "Prescribe anticoagulase medication to prevent thrombotic events."
    }

    print("Step 1: Identify the primary, most severe risk after heart valve surgery.")
    risk_explanation = "Heart valve replacement, especially with a mechanical valve, introduces a foreign surface into the circulatory system. This surface is prone to activating the body's clotting cascade, leading to the formation of blood clots (thrombi). A clot can break off and travel to the brain, causing a stroke, which is a major adverse complication."
    print(textwrap.fill(risk_explanation, width=80))
    print("-" * 80)

    print("Step 2: Evaluate each answer choice against this primary risk.")
    
    # Rationale for why other choices are less critical in this specific context
    print("\nAnalysis of other options:")
    less_critical_rationale = {
        'A & F': "Incorrect. Prophylaxis (prevention) is crucial even when asymptomatic, as the risk is high.",
        'B': "Important for comfort, but does not prevent major complications like stroke.",
        'C, D, H': "Important for long-term recovery (rehab, exercise, diet) but are not the immediate, critical step to prevent life-threatening clots.",
        'E & G': "Relate to logistics and follow-up, not the primary preventative medical intervention required upon discharge."
    }
    for key, reason in less_critical_rationale.items():
        print(f" - Options ({key}): {reason}")
    
    print("-" * 80)
    print("Step 3: Identify the action that directly addresses the primary risk.")
    
    # Explanation for the correct choice
    correct_choice_key = 'J'
    correct_choice_text = choices[correct_choice_key]
    
    conclusion = f"The only option that directly and critically addresses the risk of thrombotic events (like stroke) after heart valve surgery is the prescription of anticoagulant medication. This is the standard of care."
    print(textwrap.fill(conclusion, width=80))
    
    print("\nFinal Conclusion:")
    print(f"The most appropriate course of action is J: {correct_choice_text}")


analyze_post_op_care()