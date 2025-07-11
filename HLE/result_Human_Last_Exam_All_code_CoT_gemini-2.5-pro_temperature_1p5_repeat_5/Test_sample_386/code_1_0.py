def culture_medium_calculation():
    """
    This script demonstrates a simple calculation using the concentrations
    of key supplements mentioned in the correct cell culture protocol.
    The protocol uses Fetal Bovine Serum (FBS) and an antibiotic solution.
    """

    # Percentage of Fetal Bovine Serum (FBS) in the medium
    fbs_percentage = 10

    # Percentage of antibiotic solution in the medium
    antibiotic_percentage = 1

    # Calculate the total percentage of these two main supplements
    total_supplement_percentage = fbs_percentage + antibiotic_percentage

    print("The culture medium is composed of a basal medium plus key supplements.")
    print(f"The concentration of FBS is {fbs_percentage}%.")
    print(f"The concentration of antibiotic is {antibiotic_percentage}%.")
    print("\nCalculating the sum of these supplement percentages:")
    print(f"{fbs_percentage} + {antibiotic_percentage} = {total_supplement_percentage}")

culture_medium_calculation()