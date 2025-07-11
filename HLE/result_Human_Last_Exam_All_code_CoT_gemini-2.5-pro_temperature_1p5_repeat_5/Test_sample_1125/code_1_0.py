import random

def solve_hat_puzzle():
    """
    This function simulates the guaranteed strategy for the 12-person hat puzzle.
    """
    people = [
        "Alex", "Ben", "Charlie", "Dan", "Eva", "Frank",
        "George", "Hannah", "Ivy", "Jack", "Kyle", "Liam"
    ]
    hat_numbers = list(range(1, 13))
    random.shuffle(hat_numbers)
    
    # Each person gets a hat, but no one knows the numbers initially.
    assignments = dict(zip(people, hat_numbers))
    
    unknowns = set(people)
    knowns = {} # This will store {person: number} for those whose numbers are revealed.

    print("The team begins their strategy.")
    print("-" * 30)
    
    num_reveals = 0
    # The strategy requires 11 reveals.
    while len(unknowns) > 1:
        num_reveals += 1
        # 1. The team picks two members from the 'unknown' pool.
        p1, p2 = random.sample(list(unknowns), 2)
        
        # 2. The leader (adversarially) chooses one to reveal.
        # From the team's perspective, it doesn't matter who is chosen.
        revealed_person = random.choice([p1, p2])
        revealed_number = assignments[revealed_person]
        
        # 3. The information is now public. The revealed person moves to the 'known' group.
        knowns[revealed_person] = revealed_number
        unknowns.remove(revealed_person)
        
        print(f"Round {num_reveals}:")
        print(f"  - Pair ({p1}, {p2}) raises hands.")
        print(f"  - Leader reveals that {revealed_person}'s number is {revealed_number}.")
        print(f"  - Remaining unknowns: {len(unknowns)}")

    print("-" * 30)
    print("The information gathering is complete.")

    # 4. The 11 people whose numbers were revealed know their number.
    num_revealed = len(knowns)
    print(f"\n{num_revealed} people know their number from direct reveals:")
    for person, number in knowns.items():
        print(f"  - {person} knows their number is {number}.")

    # 5. The last unknown person deduces their number.
    last_unknown_person = list(unknowns)[0]
    all_possible_numbers = set(range(1, 13))
    known_numbers = set(knowns.values())
    last_number = (all_possible_numbers - known_numbers).pop()
    
    num_deduced = 1
    print(f"\n{num_deduced} person can deduce their number by elimination:")
    print(f"  - {last_unknown_person} deduces their number must be {last_number}, as it's the only one left.")

    # 6. Conclusion: All 12 people can determine their number.
    total_knowers = num_revealed + num_deduced
    
    print("\n--- Final Conclusion ---")
    print(f"Number of people with revealed hats: {num_revealed}")
    print(f"Number of people who can deduce their hat: {num_deduced}")
    print(f"Total number of people who can determine their number: {num_revealed} + {num_deduced} = {total_knowers}")
    
# Execute the simulation
solve_hat_puzzle()

<<<12>>>