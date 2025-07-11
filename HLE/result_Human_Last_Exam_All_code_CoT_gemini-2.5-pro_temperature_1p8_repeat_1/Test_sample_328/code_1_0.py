import sys

def solve_age_estimate():
    """
    This function provides an estimate for the age of two mating insect species
    based on their life history traits.
    """

    # Species names
    species_1_name = "Tridactylophagus tartari (Blister Beetle)"
    species_2_name = "Periplaneta americana (American Cockroach)"

    # Estimated ages based on biological data
    # Tridactylophagus tartari: Short adult lifespan, mates soon after eclosion.
    age_1_value = 3
    age_1_unit = "days"

    # Periplaneta americana: Long adult lifespan, so the average age of a mating
    # individual is expected to be much higher than its age at sexual maturity.
    age_2_value = 6
    age_2_unit = "months"

    # Explanation
    print("Step 1: Estimate the age of a mating male Tridactylophagus tartari.")
    print(f"The species {species_1_name} has a relatively short adult lifespan (weeks to a few months).")
    print("Mating typically occurs soon after eclosion. A reasonable average estimate for the male's age is therefore low.")
    print(f"-> Estimated Age = {age_1_value} {age_1_unit}\n")

    print("Step 2: Estimate the age of a mating male Periplaneta americana.")
    print(f"The species {species_2_name} has a very long adult lifespan (often over a year).")
    print("Although they mature quickly, they mate throughout their long lives. A randomly observed mating male is likely to be much older on average.")
    print(f"-> Estimated Age = {age_2_value} {age_2_unit}\n")
    
    print("--- Final Answer ---")
    print("The final combined estimate is:")
    print(f"'{species_1_name} age: {age_1_value} {age_1_unit}', '{species_2_name} age: {age_2_value} {age_2_unit}'")


solve_age_estimate()
<<<C>>>