def find_cited_lecturers():
    """
    This function provides the names of the physicists whose lectures and writings
    were foundational to Albert Einstein's 1905 paper on Brownian motion.
    """
    # The year of the paper in question.
    year = 1905

    # Historical sources indicate that Einstein's understanding of statistical mechanics
    # and thermodynamics, crucial for his paper on Brownian motion, was built
    # upon the works of Ludwig Boltzmann and Max Planck.
    lecturers = "Ludwig Boltzmann and Max Planck"

    # Construct the final answer string.
    answer = (
        f"In the development of his groundbreaking {year} paper on Brownian motion, "
        "Albert Einstein was heavily influenced by and cited the foundational work and published lectures of "
        f"{lecturers}."
    )

    print(answer)

find_cited_lecturers()