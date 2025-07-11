import random

def solve_hat_puzzle():
    """
    This script simulates the optimal strategy for the 12-person hat puzzle
    to determine the maximum number of people (N) guaranteed to find their number.
    """

    # 1. Setup the puzzle
    people = [
        'Alex', 'Ben', 'Charlie', 'Dan', 'Eva', 'Frank',
        'George', 'Hannah', 'Ivy', 'Jack', 'Kyle', 'Liam'
    ]
    numbers = list(range(1, 13))
    random.shuffle(numbers)

    # Assign a unique number to each person's hat (this is secret information)
    hats = {person: number for person, number in zip(people, numbers)}
    
    print("--- Strategy Explanation ---")
    print("The team's optimal strategy is to designate one person as the 'pivot'.")
    print("This pivot will be paired with every other person in sequence.")
    print("This simulation demonstrates the strategy in the worst-case scenario,")
    print("where the leader always reveals the number of the non-pivot person.\n")

    # This dictionary will store the information revealed to the group
    known_information = {}

    # 2. Simulate the strategy
    pivot = 'Alex'
    others = [p for p in people if p != pivot]

    print(f"--- Strategy Execution (Pivot: {pivot}) ---")
    for i, other_person in enumerate(others):
        # A pair (pivot, other_person) raises their hands.
        # The leader reveals the number of the non-pivot person.
        revealed_person = other_person
        revealed_number = hats[revealed_person]
        
        print(f"Query {i+1}: ({pivot}, {revealed_person}) -> Leader reveals {revealed_person}'s number is {revealed_number}.")
        
        # This information becomes public knowledge.
        known_information[revealed_person] = revealed_number

    print("\n--- Deduction Phase ---")
    
    # 3. Deduce the pivot's number by elimination
    print(f"After 11 queries, the numbers for 11 people are known: {known_information}")
    
    print(f"\nNow, {pivot} deduces their number.")
    all_possible_numbers = set(range(1, 13))
    revealed_numbers = set(known_information.values())
    
    # The "equation" to find the pivot's number
    missing_number_set = all_possible_numbers - revealed_numbers
    pivot_deduced_number = missing_number_set.pop()

    print("The final deduction is made by subtracting the set of revealed numbers from the set of all possible numbers:")
    # Output each number in the final equation as requested
    print(f"{{{', '.join(map(str, sorted(list(all_possible_numbers))))}}} - {{{', '.join(map(str, sorted(list(revealed_numbers))))}}} = {{{pivot_deduced_number}}}")
    
    print(f"\nSo, {pivot}'s number is {pivot_deduced_number}.")
    
    # 4. Final conclusion
    num_guaranteed = len(known_information) + 1
    
    print("\n--- Conclusion ---")
    print("Result: 11 people know their number by direct revelation, and 1 person knows their number by elimination.")
    print(f"All {num_guaranteed} people are guaranteed to determine their number correctly.")
    print(f"The largest possible value of N is {num_guaranteed}.")


solve_hat_puzzle()