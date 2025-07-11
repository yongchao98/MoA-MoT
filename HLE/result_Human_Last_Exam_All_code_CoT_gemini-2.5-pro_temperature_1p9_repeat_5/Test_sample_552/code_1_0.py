def find_cited_scientist():
    """
    Identifies and prints the name of the scientist whose lectures were cited
    by Albert Einstein in his major works of 1905.
    """
    # The scientist whose lectures were a key cited influence on Einstein's 1905 work.
    scientist = "Ludwig Boltzmann"

    # The specific work cited in his 1905 doctoral thesis.
    cited_work = "Vorlesungen über Gastheorie (Lectures on Gas Theory)"

    # While Einstein's 1905 paper on Brownian motion itself has very few explicit citations,
    # his closely related doctoral thesis, also completed in 1905, provides a clear answer.
    # In his thesis, "A New Determination of Molecular Dimensions," he directly cites the work
    # of Ludwig Boltzmann, which was published as a book derived from his lectures.
    # This work was a key foundation for Einstein's revolutionary papers on statistical mechanics.

    print(f"Albert Einstein's 1905 work, particularly his doctoral thesis, cited the lectures of: {scientist}")
    print(f"The specific work cited was '{cited_work}'.")

if __name__ == "__main__":
    find_cited_scientist()