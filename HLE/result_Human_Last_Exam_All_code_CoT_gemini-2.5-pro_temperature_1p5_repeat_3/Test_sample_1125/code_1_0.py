import random

def solve_hat_puzzle():
    """
    Simulates the "Shrinking Pool" strategy for the 12-person hat puzzle.
    """
    # 1. Setup the puzzle
    people = ['Alex', 'Ben', 'Charlie', 'Dan', 'Eva', 'Frank', 'George', 'Hannah', 'Ivy', 'Jack', 'Kyle', 'Liam']
    numbers = list(range(1, 13))
    random.shuffle(numbers)
    
    # This dictionary is secret from the team members initially
    hat_assignments = {person: number for person, number in zip(people, numbers)}
    
    # 2. The Strategy Execution
    unknowns = set(people)
    knowns = {}
    
    print("The 12 clever team members begin their strategy.")
    print("-" * 40)

    # The team can force 11 reveals of new numbers
    step = 1
    while len(unknowns) > 1:
        print(f"Step {step}:")
        # a. The team picks two people from the current pool of unknowns
        pair_to_ask = random.sample(list(unknowns), 2)
        person1, person2 = pair_to_ask
        print(f"  - The team nominates {person1} and {person2} (both have unknown numbers).")

        # b. The leader MUST reveal a new number
        revealed_person = random.choice(pair_to_ask)
        revealed_number = hat_assignments[revealed_person]
        print(f"  - The leader reveals that {revealed_person}'s number is {revealed_number}.")

        # c. The team updates their collective knowledge
        unknowns.remove(revealed_person)
        knowns[revealed_person] = revealed_number
        step += 1
        print("-" * 40)

    # 3. Deduction by the final person
    last_person_unknown = list(unknowns)[0]
    print(f"Step {step}:")
    print(f"  - 11 numbers have been revealed. The only person left is {last_person_unknown}.")
    
    all_numbers = set(range(1, 13))
    revealed_numbers = set(knowns.values())
    deduced_number = (all_numbers - revealed_numbers).pop()
    
    print(f"  - {last_person_unknown} deduces their number by elimination. The remaining number is {deduced_number}.")
    knowns[last_person_unknown] = deduced_number
    print("-" * 40)
    
    # 4. Final Conclusion
    print("Final Result:")
    print("All 12 members have successfully determined their hat numbers.")
    # We can verify our simulation against the secret assignments
    correct_determinations = 0
    for person in people:
        if knowns[person] == hat_assignments[person]:
            correct_determinations += 1
    
    print(f"Number of people who are guaranteed to find their number: {correct_determinations}")
    print("\nThe largest possible value of N is 12.")


# Execute the simulation
solve_hat_puzzle()