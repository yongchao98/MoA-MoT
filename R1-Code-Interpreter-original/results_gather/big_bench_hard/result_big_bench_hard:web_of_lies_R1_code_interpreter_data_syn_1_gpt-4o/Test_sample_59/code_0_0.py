def check_statements():
    # Assume Ka tells the truth
    ka_truth = True
    
    # Fletcher says Ka tells the truth
    fletcher_truth = ka_truth
    
    # Maybelle says Fletcher lies
    maybelle_truth = not fletcher_truth
    
    # Lorine says Maybelle lies
    lorine_truth = not maybelle_truth
    
    # Crista says Lorine tells the truth
    crista_truth = lorine_truth
    
    return crista_truth

print(check_statements())