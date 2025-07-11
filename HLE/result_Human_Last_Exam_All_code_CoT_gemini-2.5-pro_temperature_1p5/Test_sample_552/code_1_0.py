def find_cited_scientist():
    """
    Identifies and prints the name of the scientist whose lectures and work
    were cited by Albert Einstein in his 1905 paper on Brownian motion.
    """
    scientist_name = "Ludwig Boltzmann"
    paper_year = 1905
    topic = "Brownian motion"

    print(f"In his {paper_year} paper on {topic}, Albert Einstein cited the work and lectures of {scientist_name}.")

if __name__ == "__main__":
    find_cited_scientist()