def solve_literary_riddle():
    """
    This script solves a literary riddle by deducing the identity of a poet based on a series of clues.
    """
    
    # Step 1: Deconstruct the riddle into its core components.
    print("Step 1: Analyzing the clues from the riddle.")
    clues = {
        "Nationality": "English",
        "Occupation": "Poet",
        "Location": "Vienna",
        "Context": "Describing wide boulevards in a Russian translation of an essay."
    }
    for key, value in clues.items():
        print(f"- {key}: {value}")
    print("-" * 30)

    # Step 2: Identify the most likely poet based on the clues.
    # W. H. Auden is the most famous English poet with a deep connection to Vienna.
    print("Step 2: Identifying the English poet connected to Vienna.")
    poet_candidate = "W. H. Auden"
    print(f"The most famous English poet who lived and died in Vienna is {poet_candidate}.")
    print("-" * 30)

    # Step 3: Identify the likely author of the essay.
    # Joseph Brodsky, a friend of Auden, frequently wrote about him in connection with Vienna.
    print("Step 3: Identifying the likely author of the essay.")
    essayist_candidate = "Joseph Brodsky"
    print(f"The essay is very likely by {essayist_candidate}, a Russian-American admirer of Auden.")
    print("-" * 30)
    
    # Step 4: Synthesizing the information to find the final answer.
    # The reference to the poet's surname is an intellectual nod within the essay.
    print("Step 4: Explaining the final connection.")
    print("Joseph Brodsky often wrote about W. H. Auden. When describing Vienna, a city")
    print("where Auden spent his last years, Brodsky would embed a reference to him.")
    print("Therefore, in a description of Vienna's wide boulevards, the surname that appears")
    print("as a literary tribute is that of Auden.")
    print("-" * 30)

    # Step 5: Presenting the final answer.
    final_surname = "Auden"
    print("The final answer is the surname of the English poet.")
    print(f"The surname = '{final_surname[0]}'{final_surname[1]}{final_surname[2]}{final_surname[3]}{final_surname[4]}'")

solve_literary_riddle()