def generate_formula_string():
    """
    This function generates and prints the Linear Logic formula F(r, z, nz)
    for the JZDEC instruction. The formula uses symbolic names for the parameters.
    """
    # Symbolic parameters for the JZDEC instruction:
    # r: the register index (1 <= r <= n)
    # z: the next instruction label if the counter is zero
    # nz: the next instruction label if the counter is non-zero
    
    # Literals used in the encoding
    S_z = "S_z"    # Represents state z
    S_nz = "S_nz"  # Represents state nz
    C_r = "C_r"    # Represents one unit of counter r
    
    # Linear Logic connectives (using standard symbols)
    # The problem uses '&' for additive conjunction ("with")
    # and '⊸' (multimap) for linear implication.
    WITH = "&"
    LOLLIPOP = "⊸"
    
    # Construct the formula F(r, z, nz).
    # It provides a choice between two branches:
    # 1. The zero-branch: Produces the new state literal S_z.
    # 2. The non-zero-branch: An implication that consumes a C_r and produces S_nz.
    formula = f"{S_z} {WITH} ({C_r} {LOLLIPOP} {S_nz})"
    
    print("The appropriate formula F(r,z,nz) is:")
    print(formula)

generate_formula_string()