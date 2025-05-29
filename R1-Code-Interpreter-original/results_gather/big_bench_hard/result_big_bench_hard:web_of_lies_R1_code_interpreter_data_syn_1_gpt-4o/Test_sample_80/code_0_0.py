def check_statements():
    # Assume Christie tells the truth
    christie_truth = True
    
    # Jamey says Christie tells the truth
    jamey_truth = christie_truth
    
    # Michael says Jamey lies
    michael_truth = not jamey_truth
    
    # Rashida says Michael tells the truth
    rashida_truth = michael_truth
    
    # Teressa says Rashida tells the truth
    teressa_truth = rashida_truth
    
    return teressa_truth

print(check_statements())