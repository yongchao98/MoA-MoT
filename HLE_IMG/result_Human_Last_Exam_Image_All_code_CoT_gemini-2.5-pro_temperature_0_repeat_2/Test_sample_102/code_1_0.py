def get_iupac_name():
    """
    This function constructs and prints the IUPAC name of the product.
    The reaction is a Methoxy-Oxy-Cope rearrangement.
    The product is a bicyclic ketone.
    The numbers in the name are determined by IUPAC nomenclature rules for bicyclic systems.
    """
    
    # The stereochemistry is determined by the transition state geometry.
    # 'exo' indicates the methyl group points away from the C9-ketone bridge.
    part1 = "(exo)-"
    
    # Locant for the methoxy group and the double bond.
    # Numbering the bridge with substituents first gives the lowest locant set {2, 4, 9}.
    num1 = 2
    part2 = "-methoxy-"
    
    # Locant for the methyl group.
    num2 = 4
    part3 = "-methylbicyclo["
    
    # Bridge lengths of the bicyclo[X.Y.Z] system.
    num3 = 3
    part4 = "."
    num4 = 3
    part5 = "."
    num5 = 1
    part6 = "]non-"
    
    # Locant for the double bond ('en').
    num6 = 2
    part7 = "-en-"
    
    # Locant for the ketone ('one').
    num7 = 9
    part8 = "-one"
    
    # Print the name, highlighting each number as requested.
    print("The IUPAC name of the product is:")
    print(f"{part1}{num1}{part2}{num2}{part3}{num3}{part4}{num4}{part5}{num5}{part6}{num6}{part7}{num7}{part8}")

get_iupac_name()