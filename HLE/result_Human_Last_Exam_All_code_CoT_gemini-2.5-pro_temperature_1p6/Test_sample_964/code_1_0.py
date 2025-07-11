import re

def identify_final_barium_salt():
    """
    This function walks through the described chemical reactions step-by-step
    to determine the final barium salt in the flask.
    """
    
    # Step 1: Initial reaction
    print("Step 1: Mixing aqueous solutions of barium chloride (BaCl₂) and silver nitrate (AgNO₃).")
    print("A double displacement reaction occurs: BaCl₂ + 2 AgNO₃ -> Ba(NO₃)₂ + 2 AgCl")
    print("This forms solid silver chloride (AgCl) and aqueous barium nitrate (Ba(NO₃)₂).")
    flask_contents_1 = ["Ba(NO₃)₂(aq)", "AgCl(s)"]
    print(f"Flask now contains: {', '.join(flask_contents_1)}")
    
    # Step 2: First drying
    print("\nStep 2: The mixture is freeze-dried, removing water.")
    flask_contents_2 = ["Ba(NO₃)₂(s)", "AgCl(s)"]
    print(f"The solid contents are: {', '.join(flask_contents_2)}")

    # Step 3: Ammonia is added
    print("\nStep 3: Ammonia (NH₃) is added, which reacts with silver chloride.")
    print("The reaction forms a soluble complex: AgCl(s) + 2 NH₃(aq) -> [Ag(NH₃)₂]Cl(aq)")
    print("The barium nitrate simply dissolves in the aqueous ammonia.")
    flask_contents_3 = ["Ba(NO₃)₂(aq)", "[Ag(NH₃)₂]Cl(aq)"]
    print(f"Flask now contains a solution of: {', '.join(flask_contents_3)}")

    # Step 4: Second drying
    print("\nStep 4: The mixture is freeze-dried again, removing ammonia and water.")
    print("The removal of ammonia reverses the complex formation: [Ag(NH₃)₂]Cl -> AgCl + 2 NH₃")
    final_contents = ["Ba(NO₃)₂(s)", "AgCl(s)"]
    print(f"The final solid contents in the flask are: {', '.join(final_contents)}")
    
    # Conclusion
    print("\n--- Final Analysis ---")
    final_barium_salt_formula = ""
    for compound in final_contents:
        if compound.startswith("Ba"):
            # Get 'Ba(NO3)2' from 'Ba(NO₃)₂(s)' and normalize it
            final_barium_salt_formula = compound.split('(')[0].replace('₃', '3') + "(" + compound.split('(')[1].split(')')[0].replace('₃', '3') + ")" + compound.split(')')[1].split('(')[0]
            final_barium_salt_formula = "Ba(NO3)2" # a simpler way
            break
            
    print(f"The barium salt remaining in the flask is Barium Nitrate.")
    print(f"Its chemical formula is {final_barium_salt_formula}.")

    # Addressing the instruction to output numbers from the formula
    print("\nAs per the instructions, here are the numbers from the final chemical formula's equation (i.e., its subscripts):")
    # Using regex to find all numbers in the formula string "Ba(NO3)2"
    numbers = re.findall(r'\d+', final_barium_salt_formula)
    for num in numbers:
      print(num)

# Run the analysis
identify_final_barium_salt()