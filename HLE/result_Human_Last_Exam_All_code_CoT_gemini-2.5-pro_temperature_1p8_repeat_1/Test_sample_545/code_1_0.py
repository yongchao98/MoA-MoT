import pandas as pd

def solve_chemistry_problem():
    """
    This script outlines the steps to identify the major product and provides its IUPAC name.
    
    The reaction proceeds in two main steps:
    1. Thermal Sulfoxide Elimination:
       The starting material, a sulfoxide with beta-hydrogens, undergoes a syn-elimination
       upon heating.
       Reactant: Ph-S(=O)-CH2-CH2-O-C(CH3)2-CH=CH2
       Intermediate Product: CH2=CH-O-C(CH3)2-CH=CH2 (an allyl vinyl ether) + Ph-SOH
    
    2. Claisen Rearrangement:
       The allyl vinyl ether intermediate, under the same thermal conditions, immediately
       undergoes a [3,3]-sigmatropic rearrangement (Claisen rearrangement).
       Intermediate: CH2=CH-O-C(CH3)2-CH=CH2
       Final Product: O=CH-CH2-C(CH3)2-CH=CH2 (a gamma,delta-unsaturated aldehyde)
    
    3. IUPAC Naming of the Final Product:
       The final aldehyde product is named based on standard IUPAC rules.
       - The principal functional group is the aldehyde (-al).
       - The longest carbon chain including the aldehyde is 5 carbons long (pent-).
       - Numbering starts from the aldehyde carbon as C1.
       - A double bond is present between C4 and C5 (pent-4-enal).
       - Two methyl groups are on C3 (3,3-dimethyl).
    """
    
    # IUPAC name of the final product
    final_product_name = "3,3-dimethylpent-4-enal"
    
    print(f"The IUPAC name of the major product is:")
    print(final_product_name)

solve_chemistry_problem()