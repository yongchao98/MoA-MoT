def solve_riddle():
    """
    This function solves the riddle based on a phonetic pun between Russian and English.
    1. The "Ring" boulevard in Vienna is "кольцо" (kol'tso) in Russian.
    2. This sounds like the surname of the English poet, Coleridge.
    """
    
    # The surname of the English poet.
    poet_surname = "Coleridge"
    
    print("The riddle describes the Vienna Ringstrasse ('Ring Boulevard').")
    print("The Russian word for 'ring' is 'кольцо' (kol'tso).")
    print("This sounds very similar to the surname of the English poet Samuel Taylor Coleridge.")
    print("\nTherefore, the surname is:")
    
    # The final 'equation' is simply the resulting name.
    # We will print each character of the name as if it were a component of an equation.
    final_answer_components = " + ".join(list(poet_surname))
    print(f"{final_answer_components} = {poet_surname}")

solve_riddle()