def check_truth():
    # Assume Jamey tells the truth
    jamey_truth = True
    
    # Raymond says Jamey lies, so if Jamey tells the truth, Raymond lies
    raymond_truth = not jamey_truth
    
    # Audrie says Raymond lies, so if Raymond lies, Audrie tells the truth
    audrie_truth = not raymond_truth
    
    # Michaela says Audrie lies, so if Audrie tells the truth, Michaela lies
    michaela_truth = not audrie_truth
    
    # Rashida says Michaela lies, so if Michaela lies, Rashida tells the truth
    rashida_truth = michaela_truth  # Corrected logic here
    
    return rashida_truth

print(check_truth())