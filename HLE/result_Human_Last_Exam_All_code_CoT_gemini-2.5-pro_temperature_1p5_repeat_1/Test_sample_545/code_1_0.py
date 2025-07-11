import textwrap

def solve_reaction_and_name_product():
    """
    This function explains the chemical reaction and determines the IUPAC name of the major product.
    """

    print("This script will determine the IUPAC name of the major product for the given reaction.\n")

    # Step 1: Decode the reactant and the reaction pathway
    print("Step 1: Analyzing the Reaction Pathway")
    print("---------------------------------------")
    reactant_name = "((2-((2-methylbut-3-en-2-yl)oxy)ethyl)sulfinyl)benzene"
    reactant_structure = "Ph-S(=O)-CH2-CH2-O-C(CH3)2-CH=CH2"
    print(f"The starting material is: {reactant_name}")
    print(f"Its structure is: {reactant_structure}\n")

    explanation1 = """
    The reaction is heated to 180 °C. The reactant is a sulfoxide with hydrogen atoms on the
    beta-carbon (Ph-S(=O)-CH2(alpha)-CH2(beta)-...). The first step is a thermal syn-elimination
    (sulfoxide pyrolysis).
    """
    print(textwrap.dedent(explanation1))

    intermediate_structure = "CH2=CH-O-C(CH3)2-CH=CH2"
    print(f"This elimination produces an intermediate allyl vinyl ether with the structure: {intermediate_structure}")
    print("It also produces phenylsulfenic acid (Ph-SOH), which is neutralized by the sodium bicarbonate.\n")

    # Step 2: Second reaction
    print("Step 2: Claisen Rearrangement")
    print("-----------------------------")
    explanation2 = """
    The allyl vinyl ether intermediate is unstable at high temperatures and immediately undergoes
    a [3,3]-sigmatropic rearrangement, also known as the Claisen rearrangement.
    """
    print(textwrap.dedent(explanation2))

    final_product_structure = "OHC-CH2-CH2-CH=C(CH3)2"
    print(f"The rearrangement converts the ether ({intermediate_structure}) into a gamma,delta-unsaturated aldehyde.")
    print(f"The final product has the structure: {final_product_structure}\n")


    # Step 3: Naming the final product
    print("Step 3: Determining the IUPAC Name")
    print("------------------------------------")
    print(f"We will now determine the IUPAC name for the final product: {final_product_structure}\n")

    print("Naming Rules Applied:")
    print("1. Identify the principal functional group: Aldehyde (-CHO). The suffix is '-al'.")
    print("   The aldehyde carbon is assigned the number '1'.")

    print("\n2. Find the longest carbon chain that includes the aldehyde carbon.")
    print("   The chain is OHC(1)-CH2(2)-CH2(3)-CH(4)=C(5)(CH3)-CH3(6).")
    print("   The chain is 6 carbons long, so the parent name is derived from 'hexane'.")

    print("\n3. Identify and number the substituents and double bonds:")
    print("   - A double bond (ene) starts at carbon '4'. This gives the stem 'hex-4-enal'.")
    print("   - Two methyl groups are located on carbon '5'. This gives the prefix '5,5-dimethyl'.")

    print("\n4. Assemble the full name from the parts: prefix + stem.")
    
    prefix = "5,5-dimethyl"
    stem = "hex-4-enal"
    final_name = f"{prefix}{stem}"
    
    print(f"   - Prefix numbers: 5,5")
    print(f"   - Stem number: 4")
    print(f"   - Assembled name: {final_name}\n")
    
    print("--- FINAL IUPAC NAME ---")
    print(final_name)
    print("------------------------")

# Run the function
solve_reaction_and_name_product()