def find_cited_lecturer():
    """
    This function provides the name of the scientist whose lectures
    were a foundational source for Einstein's 1905 paper on Brownian motion.
    """
    scientist = "Ludwig Boltzmann"
    lecture_work = "Lectures on Gas Theory"
    paper_year = 1905
    paper_topic = "Brownian motion"

    print(f"In his {paper_year} paper on {paper_topic}, Albert Einstein's work was fundamentally based on the principles developed by {scientist}.")
    print(f"While not a direct citation in the modern sense, the paper relies heavily on the ideas from {scientist}'s influential work, famously presented in his '{lecture_work}'.")

# Execute the function to print the answer.
find_cited_lecturer()