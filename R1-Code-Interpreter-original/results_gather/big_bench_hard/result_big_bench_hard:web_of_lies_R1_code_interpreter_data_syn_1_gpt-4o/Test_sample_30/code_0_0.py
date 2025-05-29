def check_statements():
    # Assume Kandi tells the truth
    kandi_truth = True
    
    # Michaela says Kandi lies
    michaela_truth = not kandi_truth
    
    # Vernell says Michaela lies
    vernell_truth = not michaela_truth
    
    # Antwan says Vernell tells the truth
    antwan_truth = vernell_truth
    
    # Sal says Antwan lies
    sal_truth = not antwan_truth
    
    return sal_truth

print(check_statements())