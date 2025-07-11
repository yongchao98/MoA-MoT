def get_bee_common_name_info():
    """
    This function provides the common name and the year it was coined for the bee genus 'Perdita'.
    """
    
    # The genus in question is "Perdita".
    # Based on entomological sources and publications, its common name is "Fairy bees".
    common_name = "Fairy bees"
    
    # The term was used in the publication "Native bees are a rich resource in urban California gardens"
    # by Frankie, G. W., et al. in the journal California Agriculture.
    # We will use the year of this publication as the "coined" date.
    year = 2005

    # The final answer is formatted as "<name>, <year>".
    # The print statement below will construct and display this string.
    print(f"{common_name}, {year}")

get_bee_common_name_info()