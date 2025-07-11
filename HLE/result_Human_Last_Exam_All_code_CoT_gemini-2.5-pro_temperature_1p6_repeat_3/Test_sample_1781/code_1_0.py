def solve_history_puzzle():
    """
    This function analyzes the historical question and presents the correct answer with an explanation.
    """
    # Key historical facts
    monarch_name = "Philip II Augustus"
    event_description = "The change in royal title from 'King of the Franks' to 'King of France'"
    year_of_event_start = 1190
    biographer_name = "Rigord"
    monarch_epithet = "Augustus"
    year_in_correct_option = 1223 # Year of Philip II's death

    # Step-by-step explanation
    print("Step 1: Identify the monarch and the historical event.")
    print(f"The monarch in question is {monarch_name}.")
    print(f"The event is {event_description}, signifying a shift from personal to territorial rule.")
    print(f"This stylistic change began around the year {year_of_event_start}.")
    print("-" * 20)

    print("Step 2: Identify the biographer and the source of the epithet.")
    print(f"The monarch's epithet was '{monarch_epithet}'.")
    print(f"The mainstay biography and source for this epithet was written by the chronicler {biographer_name}.")
    print("-" * 20)
    
    print("Step 3: Evaluate the options and select the best fit.")
    print("We are looking for the option that combines a relevant year with the correct biographer.")
    print(f"Option C provides the pair: Year {year_in_correct_option}, Author {biographer_name}.")
    print(f"-> {biographer_name} is the correct author.")
    print(f"-> {year_in_correct_option} is the year of {monarch_name}'s death, concluding the reign during which the change occurred.")
    print("This makes it the best and most accurate choice among the options.")
    print("-" * 20)

    print("Final Answer Details:")
    print(f"Year from chosen option: {year_in_correct_option}")
    print(f"Author from chosen option: {biographer_name}")

solve_history_puzzle()