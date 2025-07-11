def find_town_name():
    """
    This function identifies the town based on visual evidence and given constraints.
    - The house is located at 435 Co Rd 1300 E, Sullivan, IL.
    - The town of Sullivan, IL has a population of ~4,326, which is less than 5,000.
    - The nearby city of Mattoon, IL has a population over 5,000 (~16,870).
    - Given the proximity and the population constraint, Mattoon is the correct answer.
    """
    town_name = "Mattoon"
    print(town_name)

find_town_name()