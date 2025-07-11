def demonstrate_jtb_regress(belief_name, justification_level=1):
    """
    This function demonstrates the infinite regress problem of JTB justification.
    To know a belief, its justification must also be known, which requires another
    justification, leading to an endless chain.
    """
    try:
        # Define the current justification in the chain, e.g., J1, J2, etc.
        current_justification_name = f"J{justification_level}"

        # Explain the current step in the regress
        print(f"Step {justification_level}: To KNOW '{belief_name}', our belief must be justified by '{current_justification_name}'.")
        print(f"   - But for '{current_justification_name}' to be valid, we must first KNOW that it is true.")
        print("-" * 20)

        # The regress: We now need to find the justification for our justification.
        # This is represented by calling the function again for the justification itself.
        demonstrate_jtb_regress(current_justification_name, justification_level + 1)

    except RecursionError:
        # Python has a recursion limit, which we use here to represent an infinite process
        print(f"Step {justification_level}: To KNOW 'J{justification_level-1}', we need 'J{justification_level}'...")
        print("\nAnd so on, ad infinitum.")
        print("This endless chain is the 'Infinite Regress Problem'.")
        print("We've hit Python's recursion limit, which simulates our inability to ever find a final, foundational justification.")

# Let's start with an initial proposition "P" that we want to know.
initial_belief = "The sun will rise tomorrow"
print("Demonstrating the 'Infinite Regress Problem' of Justification:\n")
demonstrate_jtb_regress(initial_belief)