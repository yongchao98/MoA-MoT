def check_statements():
    # Assume Elanor tells the truth
    elanor_truth = True
    
    # Ka says Elanor lies
    ka_truth = not elanor_truth
    
    # Delbert says Ka tells the truth
    delbert_truth = ka_truth
    
    # Michaela says Delbert lies
    michaela_truth = not delbert_truth
    
    # Sherrie says Michaela lies
    sherrie_truth = not michaela_truth
    
    return sherrie_truth

print(check_statements())