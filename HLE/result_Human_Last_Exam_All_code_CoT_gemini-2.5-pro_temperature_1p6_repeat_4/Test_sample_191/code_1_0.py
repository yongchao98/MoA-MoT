def explain_amyloid_engulfment():
    """
    This script explains the protein whose breakdown facilitates macrophage
    engulfment of amyloid plaques.
    """
    
    # Define the key components of the biological process
    protein_full_name = "Complement component 3"
    protein_short_name = "C3"
    active_fragment = "C3b"
    cell_type = "Macrophage / Microglia"
    target = "Amyloid"

    print(f"The protein that, when broken down, allows for the engulfment of {target} by a {cell_type} is {protein_full_name} ({protein_short_name}).")
    print("\nHere is a conceptual equation of the process:")
    
    # Print out each part of the "equation"
    print(f"Step 1: The protein '{protein_short_name}' is cleaved.")
    print(f"Step 2: This breakdown produces an active fragment, '{active_fragment}'.")
    print(f"Step 3: The '{active_fragment}' fragment coats the {target}, flagging it for removal.")
    print(f"Step 4: A '{cell_type}' recognizes '{active_fragment}' and engulfs the {target}.")

explain_amyloid_engulfment()