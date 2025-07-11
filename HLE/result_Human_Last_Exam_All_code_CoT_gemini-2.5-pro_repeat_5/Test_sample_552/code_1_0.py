def find_cited_lecturer():
    """
    This function provides the name of the person whose lectures
    Albert Einstein cited in his 1905 paper on Brownian motion.
    """
    person_name = "Ludwig Boltzmann"
    lecture_series_title = "Lectures on Gas Theory"
    
    print(f"In his 1905 paper on Brownian motion, Albert Einstein cited the lectures of {person_name}.")
    print(f"The specific work cited was his '{lecture_series_title}'.")

if __name__ == "__main__":
    find_cited_lecturer()