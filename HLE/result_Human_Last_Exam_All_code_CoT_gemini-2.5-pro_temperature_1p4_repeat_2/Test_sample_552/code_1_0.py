def find_cited_scientist():
    """
    This function provides information about the scientist whose work
    was foundational to Einstein's 1905 paper on Brownian motion.
    """
    # The scientist whose work on statistical mechanics was foundational
    # for Einstein's paper.
    scientist_name = "Ludwig Boltzmann"

    # Explanation of the context.
    explanation = (
        "In his 1905 paper on Brownian motion, 'On the motion of small particles\n"
        "suspended in a stationary liquid, as required by the molecular-kinetic\n"
        "theory of heat,' Albert Einstein built his theoretical treatment upon\n"
        "the principles of statistical mechanics.\n\n"
        "The primary architect of this field, whose work (disseminated through\n"
        "publications and lectures) was foundational for Einstein's own, was:"
    )

    # Print the information for the user.
    print(explanation)
    print("\n" + scientist_name)

if __name__ == "__main__":
    find_cited_scientist()