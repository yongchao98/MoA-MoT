import sys

# Disable writing of .pyc files for a clean execution environment
sys.dont_write_bytecode = True

def solve_braveheart_expression():
    """
    Analyzes and explains the expression pattern of the Braveheart gene
    during cardiac differentiation.
    """

    # --- Biological Information ---
    # 1. In pluripotent Embryonic Stem Cells (ESCs), Braveheart (Bvht) expression is off/low.
    # 2. During cardiac differentiation, its expression is strongly induced (INCREASES) to guide the process.
    # 3. In later stages, as cardiomyocytes mature, its expression subsides from its peak (DECREASES).

    # --- Match to Answer Choices ---
    # The answer choices describe the change over the differentiation timeline.
    # Choice C: "It increases in embryonic stem cells and decreases in differentiating heart cells."
    # - "increases in embryonic stem cells" refers to the upregulation from the ESC state.
    # - "decreases in differentiating heart cells" refers to the subsidence from the expression peak in later stages.

    chosen_answer = {
        "choice": "C",
        "text": "It increases in embryonic stem cells and decreases in differentiating heart cells."
    }

    # --- Print the step-by-step explanation ---
    print("Step 1: Analyze the change from Embryonic Stem Cells (ESCs)")
    print("As ESCs begin to differentiate into heart cells, Braveheart gene expression significantly INCREASES.")
    print("This is because its role is to initiate and control the commitment to the cardiac lineage.\n")

    print("Step 2: Analyze the expression in later-stage differentiating cells")
    print("After peaking during the differentiation process, Braveheart expression then DECREASES as the cells mature.\n")

    print("--- Conclusion ---")
    print("The complete expression dynamic involves an initial increase followed by a decrease from the peak.")
    print(f"This matches the description provided in choice '{chosen_answer['choice']}'.\n")

    # --- Final Answer Output ---
    print(f"The correct statement is: {chosen_answer['text']}")


if __name__ == "__main__":
    solve_braveheart_expression()