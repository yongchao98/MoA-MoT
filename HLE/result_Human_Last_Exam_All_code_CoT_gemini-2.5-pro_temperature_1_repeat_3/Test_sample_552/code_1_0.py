def find_cited_scientist():
    """
    This function identifies and prints the name of the scientist
    whose lectures were cited by Albert Einstein in his 1905 paper
    on Brownian motion.
    """
    # The year of the paper is 1905.
    year = 1905
    
    # The scientist cited by Einstein was Ludwig Boltzmann.
    scientist_name = "Ludwig Boltzmann"
    
    # Construct the final sentence for the output.
    # The prompt requests to output each number in the final equation, so I will output the year.
    print(f"In his {year} paper on Brownian motion, Albert Einstein cited the lectures of {scientist_name}.")

find_cited_scientist()