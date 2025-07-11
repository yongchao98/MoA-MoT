def explain_covalency():
    """
    This function explains the relationship between orbital overlap and covalency
    and applies it to the specific case of CeF6(2-) and CeCl6(2-).
    """
    
    # Statement from the problem
    overlap_in_CeF6 = "greater"
    overlap_in_CeCl6 = "lesser" # inferred
    
    # The fundamental principle
    principle = "In chemical bonding, greater orbital overlap leads to a more effective sharing of electrons between atoms, which by definition means a stronger covalent character (stronger covalency)."
    
    # Applying the principle to the given compounds
    conclusion = f"""
    Given that the orbital overlap in CeF6(2-) is {overlap_in_CeF6} than in CeCl6(2-),
    it follows from the principle of chemical bonding that CeF6(2-) will display stronger covalency.
    """
    
    print("### Principle ###")
    print(principle)
    print("\n### Application ###")
    print(conclusion)

explain_covalency()