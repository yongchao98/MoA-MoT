import fractions

def solve_pearl_riddle():
    """
    This function solves the pearl necklace riddle step-by-step.
    """
    
    # Part 1: How many pearls were there altogether?
    
    print("--- Part 1: Finding the total number of pearls ---")
    
    # The number of pearls that remained on the string
    remained_on_string = (11 * 11) - 7
    
    # Fractions of the total pearls (x) that fell
    f_floor = fractions.Fraction(1, 6)
    f_bed = fractions.Fraction(1, 5)
    f_woman = fractions.Fraction(1, 3)
    f_lover = fractions.Fraction(1, 10)
    
    # Sum of the fractions of fallen pearls
    total_fallen_fraction = f_floor + f_bed + f_woman + f_lover
    
    # The fraction of pearls that remained on the string is 1 - the fraction that fell
    remained_fraction = 1 - total_fallen_fraction
    
    # The equation is: total_pearls * remained_fraction = remained_on_string
    # So, total_pearls = remained_on_string / remained_fraction
    total_pearls = remained_on_string / remained_fraction
    
    # We want to display the full equation as requested
    print("Let 'x' be the total number of pearls on the necklace.")
    print("The equation derived from the riddle is:")
    print(f"x = ({f_floor.numerator}/{f_floor.denominator})x + ({f_bed.numerator}/{f_bed.denominator})x + ({f_woman.numerator}/{f_woman.denominator})x + ({f_lover.numerator}/{f_lover.denominator})x + {remained_on_string}")
    print("\nSolving for 'x':")
    
    # Print the value of each part
    pearls_floor = f_floor * total_pearls
    pearls_bed = f_bed * total_pearls
    pearls_woman = f_woman * total_pearls
    pearls_lover = f_lover * total_pearls
    
    print(f"  - Pearls that fell to the floor (1/6): {int(pearls_floor)}")
    print(f"  - Pearls that fell upon the bed (1/5): {int(pearls_bed)}")
    print(f"  - Pearls saved by the woman (1/3): {int(pearls_woman)}")
    print(f"  - Pearls caught by her lover (1/10): {int(pearls_lover)}")
    print(f"  - Pearls that remained upon the string: {remained_on_string}")
    
    # Summing up the parts to verify
    total_calculated = pearls_floor + pearls_bed + pearls_woman + pearls_lover + remained_on_string
    print(f"\nVerification: {int(pearls_floor)} + {int(pearls_bed)} + {int(pearls_woman)} + {int(pearls_lover)} + {remained_on_string} = {int(total_calculated)}")
    
    print(f"\nThere were {int(total_pearls)} pearls altogether on the necklace.")
    
    # Part 2: How many more are needed?
    
    print("\n--- Part 2: Calculating how many more pearls are needed ---")
    
    # Calculate total fallen pearls
    fallen_pearls = total_pearls - remained_on_string
    print(f"Total number of fallen pearls: {int(fallen_pearls)}")
    
    # Calculate how many of the fallen were found
    found_pearls = fallen_pearls * fractions.Fraction(1, 3)
    print(f"Number of fallen pearls they find back (1/3): {int(found_pearls)}")
    
    # Calculate how many are still missing/needed
    needed_pearls = fallen_pearls - found_pearls
    print(f"The number of pearls they still need is the total fallen minus the found ones ({int(fallen_pearls)} - {int(found_pearls)}).")
    
    print(f"\nThey are now gonna need {int(needed_pearls)} more pearls.")

# Run the solver
solve_pearl_riddle()