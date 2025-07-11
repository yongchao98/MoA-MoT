def get_iupac_name():
    """
    This function outlines the reaction sequence and prints the IUPAC name of the major product.
    """
    print("The reaction proceeds in two sequential steps:")
    print("\nStep 1: Sulfoxide Pyrolysis (syn-Elimination)")
    print("The starting material, ((2-((2-methylbut-3-en-2-yl)oxy)ethyl)sulfinyl)benzene, is heated.")
    print("It undergoes a thermal syn-elimination to produce an alkene and benzenesulfenic acid.")
    print("The intermediate alkene formed is an allyl vinyl ether derivative: CH2=CH-O-C(CH3)2-CH=CH2")

    print("\nStep 2: Claisen Rearrangement")
    print("The allyl vinyl ether intermediate is unstable at 180 °C and immediately undergoes a [3,3]-sigmatropic rearrangement (Claisen rearrangement).")
    print("This rearrangement is thermodynamically favorable and yields a gamma,delta-unsaturated aldehyde.")
    print("The final major product has the structure: (CH3)2C=CH-CH2-CH2-CHO")

    final_name = "5,5-dimethylhex-4-enal"
    print("\n--- Final IUPAC Name ---")
    print(f"The IUPAC name of the major product is: {final_name}")
    
    # As requested, outputting each number in the final name's equation/structure
    print("\nThe numbers used in the IUPAC name are derived from carbon positions:")
    print("Position of the two methyl groups: 5, 5")
    print("Position of the double bond (ene): 4")

if __name__ == '__main__':
    get_iupac_name()