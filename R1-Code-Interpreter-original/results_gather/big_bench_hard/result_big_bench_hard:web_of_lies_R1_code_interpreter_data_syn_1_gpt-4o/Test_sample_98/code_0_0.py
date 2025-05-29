def check_statements():
    # Assume Millicent tells the truth
    millicent_truth = True
    
    # Alexis says Millicent lies
    alexis_truth = not millicent_truth
    
    # Jaymie says Alexis lies
    jaymie_truth = not alexis_truth
    
    # Willian says Jaymie tells the truth
    willian_truth = jaymie_truth
    
    # Leda says Willian lies
    leda_truth = not willian_truth
    
    return leda_truth

print(check_statements())