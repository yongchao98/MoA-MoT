def solve_riddle():
    """
    This script solves the historical riddle by outlining the clues and connecting them to the final answer.
    """
    
    # Clues from the riddle
    poet = "Andrei Voznesensky"
    description = "a joy-discovering sailor"
    grave_mystery = "Grave was lost, but its location rediscovered in the late 1980s"
    
    # Reasoning
    print("Analyzing the clues to find the man's last name:")
    print(f"1. The poet, {poet}, wrote the libretto for the famous rock opera 'Juno and Avos'.")
    print(f"2. The opera is about the historical figure Nikolai Petrovich Rezanov (1764-1807), a Russian explorer whose journey fits the description of a '{description}'.")
    print(f"3. Rezanov died in Krasnoyarsk. His grave was destroyed and lost for decades.")
    print(f"4. The clue '{grave_mystery}' is the key. In the late 1980s, fueled by the opera's success, Rezanov's burial site was rediscovered.")
    print("\nBased on these facts, the man's last name is:")
    
    # The final answer
    last_name = "Rezanov"
    print(last_name)

solve_riddle()