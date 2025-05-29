def check_truth():
    # Assume Vernell tells the truth
    vernell_truth = True
    
    # Kandi says Vernell tells the truth
    kandi_truth = vernell_truth
    
    # Alexis says Kandi tells the truth
    alexis_truth = kandi_truth
    
    # Sima says Alexis lies
    sima_truth = not alexis_truth
    
    # Phoebe says Sima tells the truth
    phoebe_truth = sima_truth
    
    return phoebe_truth

print(check_truth())