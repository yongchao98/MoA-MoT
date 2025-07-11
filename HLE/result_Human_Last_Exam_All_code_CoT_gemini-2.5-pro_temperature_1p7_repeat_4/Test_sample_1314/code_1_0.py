def solve_riddle():
    """
    This function provides the answer to the user's question about the doppelgänger syndrome.
    """
    primary_disease = "Paranoid Schizophrenia"
    explanation = (
        "The doppelgänger syndrome is a delusional misidentification syndrome. "
        "With neuroscientific development and improvements in imaging diagnosis, a specific "
        "condition identified as the primary disease underlying this condition is {}.\n\n"
        "This disorder often involves complex delusions and hallucinations, "
        "including the experience of seeing one's own double (heautoscopy)."
    ).format(primary_disease)
    
    print(explanation)

solve_riddle()