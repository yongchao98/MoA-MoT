import random

def solve_hat_puzzle():
    """
    This script explains and demonstrates the strategy that guarantees
    all 12 team members can determine their hat number.
    """

    # 1. Setup the puzzle with 12 people and 12 unique numbers.
    people = ["Alex", "Ben", "Charlie", "Dan", "Eva", "Frank",
              "George", "Hannah", "Ivy", "Jack", "Kyle", "Liam"]
    all_possible_numbers = set(range(1, 13))

    # Assign a random, unique number from 1 to 12 to each person.
    # This represents the secret distribution of numbers.
    hat_numbers = list(range(1, 13))
    random.shuffle(hat_numbers)
    assignments = {person: number for person, number in zip(people, hat_numbers)}

    # We will analyze the strategy from one person's perspective.
    focus_person = "Alex"
    alex_actual_number = assignments[focus_person]
    other_people = [p for p in people if p != focus_person]

    print(f"The strategy: Every possible pair of team members will raise their hands.")
    print(f"Let's demonstrate the guarantee from {focus_person}'s perspective.")
    print(f"(The actual secret number on {focus_person}'s hat is {alex_actual_number}).")
    print("-" * 50)

    # 2. Simulate the worst-case reveals for the focus_person.
    # In the worst case, the leader tries to hide the focus_person's number.
    # This means revealing the other person's number in all their pairings.
    
    alex_knowledge = {} # Stores the numbers Alex learns {person: number}
    
    print(f"Simulating the 11 pairings involving {focus_person}:")
    for partner in other_people:
        # In the pair (focus_person, partner), the leader reveals the partner's number.
        revealed_person = partner
        revealed_number = assignments[revealed_person]

        # The focus_person learns this new information.
        alex_knowledge[revealed_person] = revealed_number
        print(f"- Pairing ({focus_person}, {partner}): Leader reveals {partner}'s number is {revealed_number}.")
    
    print("-" * 50)
    
    # 3. The Deduction Process.
    # Since the leader never revealed Alex's number, Alex must use deduction.
    print(f"{focus_person}'s Reasoning:")
    print("\"My number was never directly revealed.\"")
    print("\"However, I paired with all 11 other members, and in each case, their number was announced.\"")
    
    known_numbers_of_others = set(alex_knowledge.values())
    
    print(f"\"I now know the numbers of the other 11 people: {sorted(list(known_numbers_of_others))}\"")
    
    # The final equation by elimination
    deduced_set = all_possible_numbers - known_numbers_of_others
    deduced_number = deduced_set.pop()
    
    print("\"To find my number, I subtract the set of numbers I know from the total set of numbers.\"")
    print("\nFinal Equation:")
    print(f"My Number = {{ {', '.join(map(str, sorted(list(all_possible_numbers))))} }} - {{ {', '.join(map(str, sorted(list(known_numbers_of_others))))} }}")
    print(f"\n\"Therefore, my number must be {deduced_number}.\"")

    # 4. Verification
    print("-" * 50)
    if deduced_number == alex_actual_number:
        print(f"The deduction is correct! {focus_person} successfully found their number.")
    else:
        print("There was an error in the logic.")

    print("\nThis logic works for any person. Thus, all 12 are guaranteed to know.")


# Run the demonstration
solve_hat_puzzle()
