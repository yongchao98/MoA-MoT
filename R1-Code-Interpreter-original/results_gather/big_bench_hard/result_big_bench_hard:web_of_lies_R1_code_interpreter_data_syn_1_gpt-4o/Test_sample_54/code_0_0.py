def check_statements():
    # Assume Fletcher tells the truth
    fletcher_truth = True
    
    # Ka says Fletcher tells the truth
    ka_truth = fletcher_truth
    
    # Elanor says Ka lies
    elanor_truth = not ka_truth
    
    # Vernell says Elanor lies
    vernell_truth = not elanor_truth
    
    # Delfina says Vernell tells the truth
    delfina_truth = vernell_truth
    
    return delfina_truth

print(check_statements())