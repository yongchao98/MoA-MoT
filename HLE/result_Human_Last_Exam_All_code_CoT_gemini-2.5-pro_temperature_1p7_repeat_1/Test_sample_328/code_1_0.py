import sys
# This script is for illustrative purposes to explain the reasoning.
# The core task is based on biological knowledge, not computation.

def find_mating_age_estimate():
    """
    Analyzes the life cycles of Tridactylophagus tartari and Periplaneta americana
    to estimate the average age of males at mating.
    """
    
    # Species 1: Tridactylophagus tartari (Order: Strepsiptera)
    # Males of this order have a very short adult lifespan, often less than 24 hours.
    # Their sole purpose after eclosion (emerging as an adult) is to mate. They do not feed.
    # Therefore, they must mate very soon after emerging.
    tridactylophagus_age = "16 hours"

    # Species 2: Periplaneta americana (American Cockroach)
    # This species has a much longer lifespan. After its final molt into an adult,
    # it takes about 1-2 weeks to reach sexual maturity.
    # "Relatively shortly after eclosion" in the context of an adult life that can last a year
    # makes an estimate of about a month very reasonable.
    periplaneta_age = "one month"
    
    print("Based on the life cycles of the two species:")
    print("1. Tridactylophagus tartari: This is a strepsipteran, whose males have an extremely short adult lifespan (typically hours to less than a day) focused entirely on reproduction.")
    print("2. Periplaneta americana: This is the American cockroach, which takes 1-2 weeks after becoming an adult to mature sexually and lives for about a year.")
    print("\nThe best average estimate for the male's age since eclosion for each species is therefore:")
    
    # The final step is to output the two values from the best-fitting answer choice.
    print(f"\nEstimated age for Tridactylophagus tartari: {tridactylophagus_age}")
    print(f"Estimated age for Periplaneta americana: {periplaneta_age}")
    
    # The corresponding answer choice is I.
    
if __name__ == '__main__':
    # It's important to use a function to encapsulate the logic
    # and call it from the main block.
    find_mating_age_estimate()
