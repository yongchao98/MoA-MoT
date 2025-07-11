def get_histopathological_diagnosis():
    """
    Analyzes the features of the provided image and returns the diagnosis.
    
    The image shows a kidney glomerulus with the following key features:
    1. A large, cellular, crescent-shaped structure filling Bowman's space.
    2. Compression of the glomerular tuft by this crescent.
    3. Hypercellularity within the glomerular tuft.
    
    These findings are the classic histopathological hallmarks of
    a specific type of severe glomerular disease.
    """
    
    # The primary finding is the crescent formation.
    primary_finding = "A large cellular crescent in Bowman's space"
    
    # This finding leads to a specific diagnosis.
    diagnosis = "Crescentic Glomerulonephritis"
    
    print(f"Key Finding: {primary_finding}")
    print(f"Histopathological Diagnosis: {diagnosis}")

# Execute the function to print the diagnosis
get_histopathological_diagnosis()