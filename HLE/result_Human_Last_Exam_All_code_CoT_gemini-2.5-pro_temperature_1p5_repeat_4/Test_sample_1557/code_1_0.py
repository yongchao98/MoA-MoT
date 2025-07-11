def solve_history_mcq():
    """
    Analyzes the statements about the Duke of Wellington and prints the numbers of the correct ones.
    """
    # Based on historical analysis:
    # 1. Correct - Indian logistical experience was adapted for the Peninsula.
    # 2. Incorrect - Intelligence system was also used in Europe.
    # 3. Incorrect - No such standardization occurred in 1815.
    # 4. Incorrect - Factual errors about Sandhurst's founding date and Wellington's role.
    # 5. Incorrect - Logistical systems were highly influential.
    # 6. Correct - Use of local auxiliaries became standard British practice.
    # 7. Incorrect - Indian experience was crucial to Peninsular commissariat reforms.
    # 8. Correct - Flying columns tactic was developed in India and adapted elsewhere.
    # 9. Incorrect - 1813 Charter Act was driven by politics/economics, not military principles.
    # 10. Incorrect - Indian civilian administration experience was vital in the Peninsula.

    correct_options = [1, 6, 8]
    
    # Sort the numbers in ascending order
    correct_options.sort()
    
    # Print the result as a comma-separated string
    print(','.join(map(str, correct_options)))

solve_history_mcq()