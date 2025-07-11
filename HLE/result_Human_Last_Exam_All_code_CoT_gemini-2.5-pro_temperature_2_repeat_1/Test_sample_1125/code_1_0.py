import random

def solve_hat_puzzle():
    """
    Simulates the hat puzzle strategy to show that all 12 members can find their numbers.
    """
    members = ['Alex', 'Ben', 'Charlie', 'Dan', 'Eva', 'Frank', 
               'George', 'Hannah', 'Ivy', 'Jack', 'Kyle', 'Liam']
    
    # Assign unique numbers from 1 to 12 to the hats. This is the ground truth, unknown to the members initially.
    hat_numbers = list(range(1, 13))
    random.shuffle(hat_numbers)
    true_assignments = {member: number for member, number in zip(members, hat_numbers)}

    print("The 12 members have been given hats with unique numbers from 1 to 12.\n")

    # The team's knowledge
    known_assignments = {}
    unknown_members = set(members)
    explanations = {}

    # The strategy: perform 11 queries, each time reducing the unknown pool by one.
    for i in range(11):
        print(f"--- Query {i+1} ---")
        
        # Pick any two members from the pool of those with unknown numbers
        p1, p2 = random.sample(list(unknown_members), 2)
        print(f"{p1} and {p2} raise their hands.")

        # The leader reveals the number of one of them. The choice doesn't affect the outcome.
        revealed_person = random.choice([p1, p2])
        revealed_number = true_assignments[revealed_person]

        print(f"The leader reveals that {revealed_person}'s number is {revealed_number}.\n")
        
        # The information becomes common knowledge
        known_assignments[revealed_person] = revealed_number
        explanations[revealed_person] = f"My number was directly revealed by the leader as {revealed_number}."
        unknown_members.remove(revealed_person)

    # After 11 queries, only one member is left.
    last_person = unknown_members.pop()
    print("--- Final Deduction ---")
    print(f"11 members have had their numbers revealed. {last_person} is the only one left.")

    # The last person deduces their number by elimination.
    all_possible_numbers = set(range(1, 13))
    revealed_numbers = set(known_assignments.values())
    deduced_number_set = all_possible_numbers - revealed_numbers
    deduced_number = deduced_number_set.pop()

    known_assignments[last_person] = deduced_number
    explanation_text = f"I deduced my number is {deduced_number} because it is the only number from 1-12 not taken by the other 11 members."
    explanations[last_person] = explanation_text
    print(f"{last_person} deduces their number: The revealed numbers are {sorted(list(revealed_numbers))}, so my number must be {deduced_number}.\n")

    # Final result
    print("========================================")
    print("Final Result: All 12 members determined their numbers.")
    print("The final 'equation' is the complete assignment of numbers:")
    print("========================================")
    
    # Sort by member name for consistent output
    sorted_members = sorted(known_assignments.keys())
    for member in sorted_members:
        print(f"{member}: {known_assignments[member]}")

solve_hat_puzzle()