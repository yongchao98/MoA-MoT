def check_truth():
    # Assume Phoebe tells the truth
    phoebe_truth = True
    
    # Christie says Phoebe tells the truth
    christie_truth = phoebe_truth
    
    # Fletcher says Christie tells the truth
    fletcher_truth = christie_truth
    
    # Amberly says Fletcher lies
    amberly_truth = not fletcher_truth
    
    # Raymond says Amberly tells the truth
    raymond_truth = amberly_truth
    
    return raymond_truth

print(check_truth())