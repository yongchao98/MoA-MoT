def get_protein_info():
    """
    This function provides information about the protein whose inhibition
    allows for macrophage engulfment of amyloid plaques.
    """
    protein_name = "CD33"
    
    explanation = (
        "The protein is CD33 (Siglec-3). It is an inhibitory receptor on microglia "
        "(the brain's version of macrophages). When its function is 'broken down' or inhibited, "
        "it removes a brake on the microglia's phagocytic (engulfing) activity, "
        "thereby allowing for enhanced clearance of amyloid plaques."
    )
    
    print(f"Protein: {protein_name}")
    print(f"Function: {explanation}")

get_protein_info()