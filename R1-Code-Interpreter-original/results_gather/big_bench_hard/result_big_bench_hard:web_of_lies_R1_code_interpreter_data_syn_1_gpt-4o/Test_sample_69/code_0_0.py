def check_statements():
    # Assume Phoebe tells the truth
    phoebe_truth = True
    
    # Jamey says Phoebe lies
    jamey_truth = not phoebe_truth
    
    # Audrie says Jamey tells the truth
    audrie_truth = jamey_truth
    
    # Michaela says Audrie lies
    michaela_truth = not audrie_truth
    
    # Vina says Michaela lies
    vina_truth = not michaela_truth
    
    return vina_truth

print(check_statements())