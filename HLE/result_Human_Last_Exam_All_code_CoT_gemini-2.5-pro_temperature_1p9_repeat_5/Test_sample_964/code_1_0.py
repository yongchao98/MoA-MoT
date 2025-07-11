# This script determines the final barium salt after a series of chemical reactions.

def solve_chemistry_problem():
    """
    Simulates the chemical reactions step-by-step to identify the final product.
    """
    # Step 1: Barium chloride (BaCl2) and silver nitrate (AgNO3) react.
    # The reaction is a double displacement: BaCl2 + 2AgNO3 -> Ba(NO3)2 + 2AgCl
    # The barium salt is converted from Barium Chloride to Barium Nitrate.
    barium_salt = "Barium Nitrate"
    barium_salt_formula = "Ba(NO3)2"
    
    # Subsequent steps involving freeze drying and the addition/removal of ammonia
    # do not alter the chemical identity of the Barium Nitrate.
    # - Ammonia reacts with Silver Chloride, but this reaction is reversed when ammonia is removed.
    # - Barium Nitrate does not react with ammonia.
    
    # The final barium salt is therefore the one formed in the first step.
    print(f"The final barium salt in the flask is: {barium_salt}")
    
    # The primary reaction that determines the final barium salt is:
    # BaCl2 + 2AgNO3 -> Ba(NO3)2 + 2AgCl
    # To fulfill the requirement of outputting each number in the equation,
    # we will print the balanced equation with its stoichiometric coefficients.
    
    print("\nThe primary chemical equation with its coefficients is:")
    
    # Coefficients for the balanced equation
    coeff_bacl2 = 1
    coeff_agno3 = 2
    coeff_bano32 = 1
    coeff_agcl = 2
    
    print(f"{coeff_bacl2} BaCl2 + {coeff_agno3} AgNO3 -> {coeff_bano32} {barium_salt_formula} + {coeff_agcl} AgCl")

# Run the function to get the answer.
solve_chemistry_problem()
