import math

def solve_hat_puzzle():
    """
    Calculates the maximum number of team members (N) who are guaranteed
    to determine their hat number correctly.

    The solution is based on a graph theory model of the problem.
    - The 12 members are vertices in a graph.
    - The team chooses pairings, which are edges. The best strategy is to form a
      single connected component with 11 edges (a tree), like a chain.
    - An adversarial leader will choose reveals to maximize ambiguity. This corresponds
      to making the set of un-revealed members as large as possible.
    - This set of un-revealed members is a "maximum independent set" of the graph.
    - The number of people guaranteed to know their number (N) is the total number
      of members minus the size of this maximum independent set.
    """
    num_members = 12

    # For the team's optimal strategy (a path graph), the size of the
    # maximum independent set is ceil(n / 2). This is the maximum number of people
    # the leader can leave uncertain.
    max_uncertain_group_size = math.ceil(num_members / 2)

    # N is the total number of members minus those who are left uncertain.
    # The members not in the uncertain group are guaranteed to know their numbers
    # because their numbers were revealed.
    guaranteed_winners = num_members - max_uncertain_group_size

    print("The strategy guarantees a certain number of winners by modeling the problem with a graph.")
    print(f"Total number of team members: {num_members}")
    print("The team arranges pairings in a chain. This forces the leader's hand.")
    print("The number of people the leader can keep in the dark is the size of the 'maximum independent set' of this chain.")
    print(f"Calculation for the maximum number of uncertain people: ceil({num_members} / 2) = {max_uncertain_group_size}")
    print("These uncertain individuals cannot deduce their number from the remaining pool of numbers.")
    print("The number of people guaranteed to know their number (N) is the total minus the uncertain group.")
    print(f"N = {num_members} - {max_uncertain_group_size}")
    print(f"N = {guaranteed_winners}")

solve_hat_puzzle()
<<<6>>>