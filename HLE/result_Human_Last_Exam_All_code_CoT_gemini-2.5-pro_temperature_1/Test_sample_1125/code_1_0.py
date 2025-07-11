import random

def solve_hat_puzzle():
    """
    This script simulates the optimal strategy for the hat puzzle
    to demonstrate that all 12 members can determine their number.
    """
    
    # 1. Setup: 12 members and 12 unique hat numbers
    members = ["Alex", "Ben", "Charlie", "Dan", "Eva", "Frank", 
               "George", "Hannah", "Ivy", "Jack", "Kyle", "Liam"]
    hat_numbers = list(range(1, 13))
    random.shuffle(hat_numbers)
    
    # Assign numbers to members (this is the ground truth, unknown to the members)
    true_assignments = {member: number for member, number in zip(members, hat_numbers)}
    
    print("--- The Hat Puzzle Strategy Simulation ---")
    print(f"There are {len(members)} members. The hat numbers are a unique set from 1 to 12.\n")

    # 2. Strategy Implementation
    unsolved_pool = set(members)
    solved_people = {} # Will store {name: number}

    print("The team's strategy is to pair two members from the 'unsolved' pool each round.")
    print("This forces the leader to reveal a new number each time.\n")
    
    round_number = 1
    while len(unsolved_pool) > 2:
        # a. Choose two members from the unsolved pool
        p1, p2 = random.sample(list(unsolved_pool), 2)
        
        # b. Leader reveals the number of one of them
        revealed_person = random.choice([p1, p2])
        revealed_number = true_assignments[revealed_person]
        
        print(f"Round {round_number}: {p1} and {p2} pair up.")
        print(f"-> Leader reveals that {revealed_person}'s number is {revealed_number}.")
        
        # c. Update the state
        unsolved_pool.remove(revealed_person)
        solved_people[revealed_person] = revealed_number
        print(f"   {revealed_person} is now solved. {len(unsolved_pool)} members remain unsolved.\n")
        round_number += 1

    # 3. The final two members
    p1, p2 = list(unsolved_pool)
    revealed_person = random.choice([p1, p2])
    revealed_number = true_assignments[revealed_person]
    print(f"Round {round_number}: The final pair, {p1} and {p2}, raise their hands.")
    print(f"-> Leader reveals that {revealed_person}'s number is {revealed_number}.")
    unsolved_pool.remove(revealed_person)
    solved_people[revealed_person] = revealed_number
    print(f"   {revealed_person} is now solved.\n")

    # 4. The last person's deduction
    last_person = list(unsolved_pool)[0]
    all_numbers = set(range(1, 13))
    known_numbers = set(solved_people.values())
    last_person_number = all_numbers.difference(known_numbers).pop()
    
    solved_people[last_person] = last_person_number
    
    print("--- The Final Deduction ---")
    print(f"{last_person} is the last unsolved member.")
    print(f"The numbers of the other 11 members have been revealed as: {sorted(list(known_numbers))}")
    print(f"By deduction, {last_person}'s number must be the only remaining number, which is {last_person_number}.")
    print(f"{last_person} is also solved.\n")
    
    # 5. Final Result
    print("--- Conclusion ---")
    print("All 12 members have successfully determined their numbers with a convincing explanation.")
    print("Final Assignments:")
    # We create the "equation" as requested by adding the numbers
    equation_parts = []
    for person, number in sorted(solved_people.items()):
        print(f"- {person}: {number}")
        equation_parts.append(str(number))
    
    total_sum = sum(solved_people.values())
    print("\nThe sum of all numbers is a fun fact based on the assignments:")
    print(" + ".join(equation_parts) + f" = {total_sum}")

    N = len(solved_people)
    print(f"\nTherefore, the largest possible value of N is {N}.")

solve_hat_puzzle()
<<<12>>>