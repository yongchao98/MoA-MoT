import re

def solve_wittig_reaction():
    """
    Determines the product of a Wittig reaction and prints the details.
    """
    # Step 1: Define reactants based on the problem description
    aldehyde_name = "Pivalaldehyde"
    aldehyde_formula = "C5H10O"
    
    ylide_name = "(2-(2-chlorophenyl)ethylidene)triphenyl-l5-phosphane"
    ylide_formula = "C26H22ClP"
    
    # Step 2: Determine the products of the Wittig reaction
    # The reaction replaces the =O of the aldehyde with the =C part of the ylide.
    # Aldehyde part: (CH3)3C-CH=
    # Ylide part: =CH-CH2-(C6H4Cl)
    # Product: (CH3)3C-CH=CH-CH2-(C6H4Cl)
    
    # Byproduct is always Triphenylphosphine oxide
    byproduct_name = "Triphenylphosphine oxide"
    byproduct_formula = "C18H15OP"
    
    # Step 3: Name the main product and determine its formula and stereochemistry
    # Longest chain is pentene, double bond at position 2.
    # Substituents: 1-(2-chlorophenyl), 4,4-dimethyl.
    # The ylide is unstabilized, which favors the Z (cis) isomer.
    product_name = "(Z)-1-(2-chlorophenyl)-4,4-dimethylpent-2-ene"
    
    # Calculate product formula: C(4+1+1+1+6) H(9+1+1+2+4) Cl(1) = C13H17Cl
    product_formula = "C13H17Cl"
    
    # Step 4: Construct and print the output
    
    print("Wittig Reaction Analysis:")
    print("-" * 25)
    print(f"Aldehyde: {aldehyde_name} ({aldehyde_formula})")
    print(f"Ylide: {ylide_name} ({ylide_formula})")
    print("-" * 25)
    
    # Construct the final balanced equation string
    # This also fulfills the requirement to output the numbers in the equation
    equation = f"{aldehyde_formula} + {ylide_formula} -> {product_formula} + {byproduct_formula}"
    
    print("Balanced Chemical Equation:")
    print(equation)
    
    # To satisfy the "output each number in the final equation" requirement, let's also print them explicitly
    numbers = re.findall(r'\d+', equation)
    print(f"The numbers in the equation are: {', '.join(numbers)}")

    print("\nProduct Information:")
    print(f"The major organic product is: {product_name}")
    print(f"Its molecular formula is: {product_formula}")


solve_wittig_reaction()