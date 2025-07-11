def calculate_risk_factor_change():
    """
    This function demonstrates the effect of a fold increase on a baseline value,
    in the context of antigen presentation probability.
    """
    # Let's represent the baseline probability of presentation as a relative value of 1.
    baseline_presentation_probability = 1
    
    # The variant increases this probability by a factor of 1000.
    fold_increase = 1000
    
    # Calculate the new probability of presentation.
    new_presentation_probability = baseline_presentation_probability * fold_increase
    
    print("This scenario describes a significant increase in the probability of presenting a disease-causing self-antigen.")
    print("If the original relative probability of presentation is {}, and it increases by a factor of {}, the new relative probability becomes:".format(baseline_presentation_probability, fold_increase))
    # Print the equation with all the numbers
    print("{} * {} = {}".format(baseline_presentation_probability, fold_increase, new_presentation_probability))
    print("\nThis massive increase in presenting the self-antigen significantly raises the likelihood of an autoimmune response, thus increasing the disease risk.")

calculate_risk_factor_change()