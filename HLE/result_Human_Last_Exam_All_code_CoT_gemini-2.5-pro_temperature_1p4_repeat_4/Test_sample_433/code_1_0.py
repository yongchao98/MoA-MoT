import textwrap

def get_ampule_bud():
    """
    Explains the Beyond-Use Date (BUD) for a single-dose ampule according to sterile compounding standards.
    """
    container_type = "single-dose ampule"
    environment = "sterile environment (e.g., ISO Class 5)"
    
    # According to USP <797>, an ampule is a critical distinction.
    # Once an ampule is broken open, it cannot be resealed. It becomes an open system.
    # Therefore, it cannot be stored.
    
    bud_explanation = f"""
    Container Type: {container_type}
    Opened In: {environment}
    ---------------------------------------------------------------------
    Rule: According to USP <797> guidelines for sterile compounding, once a {container_type} is opened, it cannot be resealed. It is considered an open system, and its contents are exposed to the surrounding air.

    Conclusion: The contents must be used immediately. There is no storage time for an opened ampule. Any unused portion must be discarded immediately after the initial use.

    Final Answer: The Beyond-Use Date (BUD) is for IMMEDIATE USE ONLY.
    """
    
    # The rule for single-dose VIALS is different (up to 12 hours in an ISO 5 environment),
    # but that does not apply here. There is no equation or number for the ampule's BUD, 
    # as it should not be stored at all.

    print(textwrap.dedent(bud_explanation).strip())

get_ampule_bud()