import random

def solve_hat_puzzle():
    """
    Simulates the optimal strategy for the 12-person hat puzzle
    and demonstrates that all 12 can determine their number.
    """
    people = ["Alex", "Ben", "Charlie", "Dan", "Eva", "Frank", "George", "Hannah", "Ivy", "Jack", "Kyle", "Liam"]
    
    # Let's create a random, hidden assignment of hats for the simulation.
    # The team does not have access to this information.
    hat_numbers = list(range(1, 13))
    random.shuffle(hat_numbers)
    secret_assignment = {person: hat for person, hat in zip(people, hat_numbers)}

    # The team's knowledge base
    known_hats = {}  # key: person, value: hat_number
    unknowns = list(people)
    query_log = []
    round_num = 1

    print("--- The Strategy Execution Begins ---")
    while len(unknowns) > 1:
        print(f"\nRound {round_num}: {len(unknowns)} members with unknown hats.")
        
        # Create pairs from the list of members with unknown hats
        temp_unknowns = list(unknowns)
        random.shuffle(temp_unknowns) # The order of pairing doesn't matter
        pairs = []
        while len(temp_unknowns) >= 2:
            p1 = temp_unknowns.pop(0)
            p2 = temp_unknowns.pop(0)
            pairs.append((p1, p2))

        # For each pair, have the leader reveal one hat number
        for p1, p2 in pairs:
            # The leader reveals one person's hat number. The team can't control which.
            # We'll randomly pick one to simulate the leader's choice.
            revealed_person = random.choice([p1, p2])
            revealed_number = secret_assignment[revealed_person]

            # All team members learn this new information
            known_hats[revealed_person] = revealed_number
            unknowns.remove(revealed_person)
            
            log_entry = f"Pair ({p1}, {p2}) -> Leader revealed that {revealed_person}'s hat has the number {revealed_number}."
            query_log.append(log_entry)
            print(log_entry)
        
        round_num += 1
    
    # The final deduction for the last person
    print(f"\nThere is only one member left with an unknown hat: {unknowns[0]}.")
    last_unknown_person = unknowns[0]
    all_hats_set = set(range(1, 13))
    revealed_hats_set = set(known_hats.values())
    deduced_hat_number = (all_hats_set - revealed_hats_set).pop()
    
    known_hats[last_unknown_person] = deduced_hat_number
    print(f"{last_unknown_person} can deduce their number by elimination.")

    print("\n--- Final Results: All 12 members know their numbers! ---")
    
    # Sort for consistent output
    sorted_people = sorted(people)
    for person in sorted_people:
        number = known_hats[person]
        # Determine the reason for each person's knowledge
        if person == last_unknown_person:
            explanation = f"My hat number must be {number}. All other 11 hat numbers were revealed, and this is the only number left in the set from 1 to 12."
        else:
            explanation = f"My hat number is {number}. This was revealed by the leader."
            
        print(f"{person}: I know my number is {number}. Reason: {explanation}")

solve_hat_puzzle()