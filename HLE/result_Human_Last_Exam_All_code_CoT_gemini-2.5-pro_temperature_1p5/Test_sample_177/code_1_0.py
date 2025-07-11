def solve_immunology_question():
    """
    This function analyzes the molecular requirements for a T cell to act as an
    Antigen-Presenting Cell (APC) and prints the correct receptor choice.
    """

    explanation = """
    T cells are normally responders in an immune reaction, activated by professional Antigen-Presenting Cells (APCs).
    APCs, such as dendritic cells, specialize in processing external threats (exogenous antigens) and presenting them to helper T cells.
    The specific receptor that APCs use for this task is the Major Histocompatibility Complex (MHC) class II molecule.

    Therefore, to engineer a T cell to gain the function of an APC, it must be equipped with the primary tool for this job.
    Expressing MHC class II would allow the T cell to present antigens it has taken up to other T cells, fulfilling the role of an APC.
    """

    answer_choice = "E"
    answer_text = "MHC class II"

    print("--- Analysis of the Cell Engineering Task ---")
    print(explanation)
    print(f"Conclusion: The correct receptor is {answer_choice}. {answer_text}")

solve_immunology_question()