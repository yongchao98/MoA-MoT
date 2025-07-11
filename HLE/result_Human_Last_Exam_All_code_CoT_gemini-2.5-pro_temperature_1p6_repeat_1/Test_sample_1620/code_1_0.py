def check_disclosure_ethics():
    """
    Analyzes whether a lawyer's disclosure of confidential information
    was ethically permissible based on the public safety exception.
    The boolean values are represented as 1 (True) or 0 (False) for the equation.
    """

    # --- Step 1: Evaluate the conditions for the public safety exception ---

    # Was there a reasonable belief of a threat? Eric's specific words ("gun", "tonight") make it reasonable.
    reasonable_belief_of_threat = True
    
    # Was the threat of death or serious bodily harm? "end things for both of us" with a gun is.
    risk_of_death_or_serious_harm = True

    # Was the risk imminent? "tonight" is imminent.
    risk_is_imminent = True
    
    # Was the disclosure necessary to prevent the harm? Yes, contacting police was the necessary step.
    disclosure_was_necessary = True
    
    # Convert booleans to integers for the "equation" output
    val1 = int(reasonable_belief_of_threat)
    val2 = int(risk_of_death_or_serious_harm)
    val3 = int(risk_is_imminent)
    val4 = int(disclosure_was_necessary)
    
    # --- Step 2: Determine if disclosure is permissible ---
    # The ethical rule requires all conditions to be met (a logical AND).
    disclosure_is_permissible = (
        reasonable_belief_of_threat and 
        risk_of_death_or_serious_harm and 
        risk_is_imminent and
        disclosure_was_necessary
    )

    result = int(disclosure_is_permissible)

    # --- Step 3: Print the analysis and conclusion ---
    print("Analyzing James's decision based on the Public Safety Exception to confidentiality.")
    print("The logical 'equation' is that all conditions must be true.")
    print(f"Condition 1 (Reasonable Belief): {val1}")
    print(f"Condition 2 (Risk of Death/Serious Harm): {val2}")
    print(f"Condition 3 (Imminent Risk): {val3}")
    print(f"Condition 4 (Disclosure was Necessary): {val4}")
    print("----------------------------------------------------------")
    print(f"Final logical check: {val1} AND {val2} AND {val3} AND {val4} results in: {result} (where 1 is True)")
    print("----------------------------------------------------------")
    
    if disclosure_is_permissible:
        print("\nConclusion: James's disclosure was ethically permissible.")
        print("This corresponds to answer choice E: Yes, because James reasonably believed that there was an imminent threat of death or serious bodily harm, and disclosure was necessary to prevent the possible harm.")
    else:
        print("\nConclusion: James's disclosure was not ethically permissible.")

check_disclosure_ethics()