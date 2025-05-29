def check_statements():
    # Assume Ryan lies
    ryan_lies = True
    
    # Kristian says Ryan lies
    kristian_truth = ryan_lies
    
    # Teressa says Kristian tells the truth
    teressa_truth = kristian_truth
    
    # Amberly says Teressa tells the truth
    amberly_truth = teressa_truth
    
    # Jamey says Amberly tells the truth
    jamey_truth = amberly_truth
    
    return jamey_truth

print(check_statements())