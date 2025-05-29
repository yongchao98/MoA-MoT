def check_statements():
    # Assume Teressa lies
    teressa_truth = False
    
    # If Teressa lies, Bernita tells the truth
    bernita_truth = not teressa_truth
    
    # If Bernita tells the truth, Ka tells the truth
    ka_truth = bernita_truth
    
    # If Ka tells the truth, Raymond lies
    raymond_truth = not ka_truth
    
    # If Raymond lies, Jamey does not lie
    jamey_truth = not raymond_truth
    
    # Check if the assumptions are consistent
    return not teressa_truth and bernita_truth and ka_truth and not raymond_truth and jamey_truth

print(check_statements())