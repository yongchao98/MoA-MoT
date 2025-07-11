def find_cited_lecturer():
    """
    Identifies the scientist whose foundational work, often disseminated
    through lectures and publications, was cited by Albert Einstein in his 
    1905 paper on Brownian motion.
    """
    # Einstein's 1905 paper on Brownian motion was a direct application of the
    # principles of statistical mechanics.
    # The primary architect of this field, whose work Einstein studied and built upon,
    # was Ludwig Boltzmann. While Einstein developed his specific theory independently
    # and didn't cite lectures in the modern sense, the intellectual foundation
    # came from Boltzmann.
    cited_scientist = "Ludwig Boltzmann"
    paper_year = 1905

    print(f"In his paper on Brownian motion in the year {paper_year}, Albert Einstein did not cite a specific set of lectures in the text.")
    print("However, the paper is fundamentally based on the principles of statistical mechanics, a field pioneered by the physicist whose work Einstein studied extensively.")
    print(f"That physicist was: {cited_scientist}")

find_cited_lecturer()