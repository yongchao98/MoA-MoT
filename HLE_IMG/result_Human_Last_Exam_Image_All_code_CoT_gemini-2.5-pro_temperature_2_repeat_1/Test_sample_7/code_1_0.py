def solve_formula():
    """
    This function derives and prints the formula for the second voltage plateau.
    
    Nomenclature:
    V1: Voltage of the first plateau (S2->S1 transition) = 0.09 V
    V2: Voltage of the second plateau (S3->S2 transition) = 0.13 V
    mu_k: Chemical potential of Li+ in (almost full) stage k
    e: elementary charge
    """
    
    # The problem asks for the formula for the second plateau, V2.
    # From the reasoning in the plan: V2 = V1 + (mu_2 - mu_3) / e
    # We are asked to output the final equation including numbers.
    
    v1_value = 0.09
    
    # The formula expresses V2. We will print the equation representing this relationship.
    # Note that the prompt asks for "the final equation", not a calculation.
    # Python symbols μ and e cannot be executed, so we print a string.
    
    print(f"The formula for the second plateau (V_2) is derived from the relationships between the plateaus:")
    print(f"V_1 = -(μ_2 - μ_Li) / e = {v1_value} V")
    print(f"V_2 = -(μ_3 - μ_Li) / e")
    print("\nBy substituting μ_Li from the first equation into the second, we get:")
    print(f"V_2 = (μ_2 - μ_3)/e + {v1_value}")

solve_formula()