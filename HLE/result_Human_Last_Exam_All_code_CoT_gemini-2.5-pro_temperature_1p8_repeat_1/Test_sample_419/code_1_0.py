import textwrap

def explain_experimental_control():
    """
    Explains the role of the anti-flag antibody in the described experiment
    and identifies the correct answer choice.
    """
    
    explanation = """
    In the described experiment, the main goal is to demonstrate that an antibody's binding to MUC1 is specifically blocked by its sugar target, GalNAc.
    
    The Problem:
    A high concentration of an inhibitor like 500 mM GalNAc can have unintended 'off-target' effects on the cells. A critical concern is that it could cause the cells to internalize the MUC1 protein from their surface. If this happened, the anti-MUC1 antibody would fail to bind simply because its target is no longer on the cell surface, not because of specific competitive inhibition. This would invalidate the experiment's conclusion.
    
    The Control:
    The anti-flag antibody provides the necessary control. It binds to the FLAG tag on the MUC1 protein, an interaction that is completely independent of the TN antigen (the sugar).
    - Its purpose is to act as a 'loading control' for surface MUC1. By measuring the signal from the anti-flag antibody, you can confirm whether the total amount of MUC1 on the cell surface remains constant between the control (PBS) and inhibitor (GalNAc) conditions.
    
    The Timing:
    The anti-flag antibody detects its target epitope directly on the MUC1 protein, just like the primary anti-MUC1-TN antibody does. Therefore, it is a primary antibody and must be added during the primary antibody incubation step.
    
    Conclusion:
    The anti-flag antibody is essential. It must be added with the primary antibodies to verify that the 500 mM GalNAc treatment has not caused a non-specific change in the surface expression of the MUC1 protein.
    
    This reasoning directly supports answer choice C.
    """
    
    print(textwrap.dedent(explanation).strip())
    
    # Final Answer
    final_answer = "<<<C>>>"
    print(f"\nFinal Answer: {final_answer}")

# Run the explanation
explain_experimental_control()