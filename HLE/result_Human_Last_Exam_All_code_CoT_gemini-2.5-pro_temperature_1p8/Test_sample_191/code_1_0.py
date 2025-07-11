def solve_biology_question():
    """
    This function provides the answer to the user's question about
    the protein that, when broken down, facilitates macrophage
    engulfment of amyloid.
    """
    protein_name = "Pentraxin 3 (PTX3)"
    explanation = (
        "In its full-length form, PTX3 binds to amyloid-beta and inhibits its uptake by microglia (the brain's macrophages).\n"
        "However, when PTX3 is 'broken down' or cleaved by enzymes, its resulting N-terminal fragment acts as a signal to promote the engulfment of amyloid by these same cells."
    )

    print(f"The protein is: {protein_name}")
    print("\nExplanation:")
    print(explanation)

solve_biology_question()