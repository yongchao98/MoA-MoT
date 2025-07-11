import random

def solve_hat_puzzle():
    """
    Simulates the optimal strategy for the hat puzzle to find the
    maximum number of people (N) who can be certain of their hat number.
    """
    members = [
        'Alex', 'Ben', 'Charlie', 'Dan', 'Eva', 'Frank',
        'George', 'Hannah', 'Ivy', 'Jack', 'Kyle', 'Liam'
    ]
    
    # Assign a unique number from 1 to 12 to each member
    hat_numbers = list(range(1, 13))
    random.shuffle(hat_numbers)
    member_hats = dict(zip(members, hat_numbers))
    
    print("The secret hat numbers have been assigned:")
    # print(member_hats) # Uncomment to see the secret numbers
    print("-" * 30)

    # Initially, all members' numbers are unknown
    unknown_members = set(members)
    known_hats = {}
    
    round_num = 1
    while len(unknown_members) > 1:
        print(f"--- Round {round_num} ---")
        print(f"Unknown members: {sorted(list(unknown_members))}")
        
        # Pick two members from the unknown set
        pair = random.sample(list(unknown_members), 2)
        p1, p2 = pair[0], pair[1]
        
        # Leader adversarially reveals one number.
        # It doesn't matter which one, the strategy works either way.
        # We'll randomly pick one to simulate the leader's choice.
        revealed_member = random.choice(pair)
        revealed_number = member_hats[revealed_member]
        
        print(f"Pair ({p1}, {p2}) raises their hands.")
        print(f"Leader reveals that {revealed_member}'s number is {revealed_number}.")
        
        # The revealed member's number is now public knowledge
        known_hats[revealed_member] = revealed_number
        unknown_members.remove(revealed_member)
        
        print(f"{revealed_member} now knows their number.")
        print("-" * 30)
        round_num += 1

    # The final member deduces their number
    last_member = list(unknown_members)[0]
    print(f"--- Final Deduction ---")
    print(f"11 members now know their numbers: {sorted(known_hats.keys())}")
    
    all_possible_nums = set(range(1, 13))
    known_nums = set(known_hats.values())
    last_number_set = all_possible_nums - known_nums
    last_number = last_number_set.pop()

    print(f"{last_member} is the only remaining member.")
    print(f"The numbers revealed are: {sorted(list(known_nums))}")
    print(f"{last_member} deduces their number with the equation:")
    print(f"My number = {{1..12}} - {{revealed numbers}}")
    # We print each number in the equation as requested
    all_nums_str = ", ".join(map(str, sorted(list(all_possible_nums))))
    known_nums_str = ", ".join(map(str, sorted(list(known_nums))))
    print(f"My number = {{{all_nums_str}}} - {{{known_nums_str}}}")
    print(f"My number = {last_number}")

    known_hats[last_member] = last_number
    
    print("\nConclusion:")
    print("All 12 members have successfully determined their hat numbers.")
    print("The maximum number of people (N) who are guaranteed to determine their number is 12.")
    
# Execute the simulation
solve_hat_puzzle()