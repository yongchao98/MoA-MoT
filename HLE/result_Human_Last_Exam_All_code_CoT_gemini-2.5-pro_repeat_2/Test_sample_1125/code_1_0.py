import random

def solve_hat_puzzle():
    """
    Simulates the hat puzzle strategy to show all 12 members can find their number.
    """
    # 1. Initialization
    names = ["Alex", "Ben", "Charlie", "Dan", "Eva", "Frank", "George", "Hannah", "Ivy", "Jack", "Kyle", "Liam"]
    all_numbers = set(range(1, 13))
    
    # Assign a unique, random number from 1-12 to each person's hat
    # This is the ground truth, unknown to the team members at the start.
    numbers_to_assign = list(all_numbers)
    random.shuffle(numbers_to_assign)
    hats = {name: numbers_to_assign.pop() for name in names}
    
    print("--- The Hat Puzzle ---")
    print(f"There are {len(names)} members. A unique number from 1 to {len(names)} is in each hat.\n")

    # The set of people whose hat numbers are not yet known.
    unknown_people = set(names)
    # The dictionary of people whose hat numbers have been revealed.
    known_people_and_numbers = {}

    # 2. Simulation of the Reveal Process
    print("--- The Reveal Process Begins ---")
    # The team will perform 11 queries to reveal 11 numbers.
    for i in range(len(names) - 1):
        # Pick two distinct people from the current set of unknowns.
        p1, p2 = random.sample(list(unknown_people), 2)
        
        # The leader chooses one person from the pair to reveal their number.
        # This choice can be random, or adversarial. The outcome is the same.
        revealed_person = random.choice([p1, p2])
        revealed_number = hats[revealed_person]
        
        print(f"Round {i+1}: The pair ({p1}, {p2}) raises their hands.")
        print(f" -> The leader reveals that {revealed_person}'s number is {revealed_number}.\n")
        
        # Update the collective knowledge
        known_people_and_numbers[revealed_person] = revealed_number
        unknown_people.remove(revealed_person)

    # 3. Final Deduction and Announcement
    print("--- The Reveal Process Ends. Announcing the Numbers ---")

    # The last person whose number was not revealed.
    last_person = unknown_people.pop()
    
    # The last person deduces their number.
    revealed_numbers = set(known_people_and_numbers.values())
    
    # Perform the deduction "equation"
    missing_number_set = all_numbers - revealed_numbers
    deduced_number = missing_number_set.pop()

    # Store the final result for the last person
    final_results = known_people_and_numbers.copy()
    final_results[last_person] = deduced_number

    # Print the reasoning for each person
    for name in sorted(names):
        number = final_results[name]
        if name in known_people_and_numbers:
            print(f"- {name}'s number is {number}. Reasoning: 'My number was directly revealed by the leader.'")
        else:
            print(f"- {name}'s number is {number}. Reasoning: 'The set of all numbers is {sorted(list(all_numbers))}. "
                  f"The set of revealed numbers for others is {sorted(list(revealed_numbers))}. "
                  f"The only missing number is {number}, so it must be mine.'")
    
    print("\nConclusion: All 12 members were able to determine their number with a convincing explanation.")


if __name__ == '__main__':
    solve_hat_puzzle()
