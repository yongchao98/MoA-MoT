def solve_hat_puzzle():
    """
    Solves the hat puzzle by reasoning from the blind person H's perspective.

    H's statement implies there are two possible scenarios for hat types that
    produce the observed K/DK sequence, and these scenarios differ only for Alice.
    """

    people = ['A', 'B', 'C', 'D', 'E', 'F', 'G']

    # The people who said "I know"
    knowers = {'A', 'C', 'E', 'G'}
    # The people who said "I don't know"
    dont_knowers = {'B', 'D', 'F'}

    # Scenario 1: 3 Number hats, 4 Color hats.
    # A, C, E get the Number hats. G knows for another reason (e.g., last to speak).
    scenario1_types = {}
    for person in people:
        if person in {'A', 'C', 'E'}:
            scenario1_types[person] = 'Number'
        else:
            scenario1_types[person] = 'Color'

    # Scenario 2: 4 Number hats, 3 Color hats.
    # All knowers (A, C, E, G) have Number hats.
    scenario2_types = {}
    for person in people:
        if person in knowers:
            scenario2_types[person] = 'Number'
        else:
            scenario2_types[person] = 'Color'

    print("Analyzing two plausible scenarios based on H's deduction:")
    print("-" * 40)
    print(f"{'Person':<10} | {'Scenario 1 Type':<18} | {'Scenario 2 Type':<18}")
    print("-" * 40)

    alice = None
    for person in people:
        type1 = scenario1_types[person]
        type2 = scenario2_types[person]
        
        print(f"{person:<10} | {type1:<18} | {type2:<18}")

        if type1 != type2:
            alice = person

    print("-" * 40)
    print(f"H can determine the hat type for everyone except for one person.")
    print(f"The person whose hat type is ambiguous between the two scenarios is {alice}.")
    print(f"Therefore, Alice is {alice}.")

solve_hat_puzzle()