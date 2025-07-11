def find_successor():
    """
    This function provides the answer to the historical question about the voivode of Pskov.
    The query is knowledge-based, not computational.
    """
    reasoning = (
        "In 1700, at the start of the Great Northern War, Pskov became a crucial military headquarters. "
        "Peter the Great appointed Field Marshal Boris Petrovich Sheremetev as the commander-in-chief of the forces based there. "
        "In this role, he assumed supreme military and civil authority in the region, effectively succeeding the previous voivode, Ivan Ivanovich Golovin."
    )
    
    answer_choice = "A"
    full_name = "Boris Petrovich Sheremetev"
    
    print("Reasoning:")
    print(reasoning)
    print("\nTherefore, the correct answer is:")
    print(f"{answer_choice}. {full_name}")

find_successor()