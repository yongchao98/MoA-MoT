def solve_t_cell_puzzle():
    """
    This function explains the reasoning for choosing the correct receptor
    to make a T cell function as an antigen-presenting cell (APC) and prints the final answer.
    """
    explanation = """
    1.  **Function of Antigen-Presenting Cells (APCs):** Professional APCs (like macrophages and dendritic cells) are specialized in taking up antigens from the environment (exogenous antigens), processing them into smaller peptides, and presenting these peptides to activate helper T (CD4+) cells.

    2.  **The Key Receptor for Presentation:** The molecule responsible for presenting these exogenous antigen peptides is the Major Histocompatibility Complex (MHC) class II receptor. This is the defining feature of a professional APC.

    3.  **T Cell Characteristics:** T cells do not normally express MHC class II. They are the cells that *recognize* antigens presented by APCs. While T cells do express MHC class I (like most nucleated cells), this is for presenting *internal* (endogenous) antigens, not for the APC's primary role of presenting *external* threats.

    4.  **Evaluating the Options:**
        *   **CD80/CD86:** These are co-stimulatory molecules that provide a second signal for T cell activation but do not present the antigen.
        *   **MHC class I:** T cells already have this for presenting internal antigens.
        *   **TIM-4:** This receptor is involved in clearing apoptotic cells but is not the primary antigen-presenting molecule.
        *   **MHC class II:** Engineering a T cell to express MHC class II would grant it the essential ability to present exogenous antigens to helper T cells, allowing it to function as an APC.

    Therefore, MHC class II is the correct answer.
    """
    print(explanation)
    final_answer = "E"
    print(f"The correct receptor is MHC class II.")

solve_t_cell_puzzle()
<<<E>>>