def check_statements():
    # Assume Teressa tells the truth
    teressa_truth = True
    
    # If Teressa tells the truth, Bernita lies
    bernita_truth = not teressa_truth
    
    # If Bernita lies, Ka does not tell the truth
    ka_truth = not bernita_truth
    
    # If Ka does not tell the truth, Raymond does not lie
    raymond_truth = not ka_truth
    
    # If Raymond does not lie, Jamey lies
    jamey_truth = not raymond_truth
    
    # Check if the assumptions are consistent
    return teressa_truth and not bernita_truth and not ka_truth and raymond_truth and not jamey_truth

print(check_statements())