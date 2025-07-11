def print_relationship():
    """
    This function prints the derived relationship between the concentration of
    3-Hydroxypropionate ([B]) and PEP ([F]).
    """
    initial_substrate = "[B]"
    final_product = "[F]"
    # Rate constants for the direct pathway from B to F
    # B -k2-> Malonyl-CoA -k3-> Acetyl-CoA -k4-> Pyruvate -k5-> F
    k2 = "k2"
    k3 = "k3"
    k4 = "k4"
    k5 = "k5"

    # The concentration of the final product [F] is proportional to the
    # initial substrate [B] and the product of the rate constants along the direct path.
    print(f"The relationship between [B] and [F] is determined by the direct pathway:")
    print(f"3-Hydroxypropionate -{k2}-> Malonyl-CoA -{k3}-> Acetyl-CoA -{k4}-> Pyruvate -{k5}-> PEP")
    print("\nResulting Proportionality Expression:")
    print(f"{final_product} ∝ {initial_substrate} * {k2} * {k3} * {k4} * {k5}")

print_relationship()