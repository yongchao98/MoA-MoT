def solve_covalency_question():
    """
    Solves the relative covalency problem based on orbital overlap.
    """
    # Premise provided in the problem statement
    # We represent the relationship as a string for logical processing.
    relative_overlap = "greater"  # CeF6(2-) has greater overlap than CeCl6(2-)

    # The fundamental chemical principle connecting orbital overlap and covalency.
    principle = "Greater orbital overlap leads to stronger covalency."

    # Logically derive the conclusion
    if relative_overlap == "greater":
        conclusion = "stronger"
    else:
        # This case is not the one described, but shows the full logic.
        conclusion = "weaker"

    # Output the step-by-step reasoning, as if it were an equation.
    print("Here is the logical deduction:")
    print("Step 1 (Premise): Orbital overlap in CeF6(2-) is " + relative_overlap + " than in CeCl6(2-).")
    print("Step 2 (Principle): " + principle)
    print("Step 3 (Conclusion): Therefore, CeF6(2-) displays " + conclusion + " covalency compared to CeCl6(2-).")

solve_covalency_question()