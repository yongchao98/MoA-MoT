def check_statements():
    # Assume Elanor lies
    elanor_lies = True
    
    # Tamika says Elanor lies
    tamika_truth = elanor_lies
    
    # Ka says Tamika tells the truth
    ka_truth = tamika_truth
    
    # Shaunda says Ka lies
    shaunda_truth = not ka_truth
    
    # Conception says Shaunda lies
    conception_truth = not shaunda_truth
    
    return conception_truth

print(check_statements())