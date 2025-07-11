def solve_graphene_puzzle():
    """
    This function determines the order of simulations based on the given conditions.
    The logic is based on the analysis of the tight-binding parameters from the provided plots.
    """

    # From visual analysis and tight-binding formulas, we deduced the approximate parameters:
    # Sim 1: t=2.5, s=1/6
    # Sim 2: t=3.0, s=0
    # Sim 3: t=2.5, s=1/6
    # Sim 4: t=2.5, s=-1/6
    # Note: Sim 1 and Sim 3 appear identical, creating ambiguity.

    # Assigning simulations to conditions based on analysis:
    # Condition 2: minimum |s| -> Sim 2 (s=0 is clearly the minimum magnitude).
    # Condition 3: unique sign(s) -> Sim 4 (s<0 is the unique sign).
    # Conditions 1 (min t) and 4 (max s) are ambiguous between Sim 1 and Sim 3.
    # To resolve the ambiguity and provide a unique mapping as requested by the problem format,
    # we make an arbitrary but necessary choice:
    # - Assign 'maximum s' to the simulation with the lower index (1).
    # - Assign 'minimum t' to the simulation with the higher index (3).

    sim_for_min_t = 3
    sim_for_min_s_mag = 2
    sim_for_unique_sign = 4
    sim_for_max_s = 1

    # The final answer is the sequence of simulation indices ordered by the conditions (1 to 4).
    result = [sim_for_min_t, sim_for_min_s_mag, sim_for_unique_sign, sim_for_max_s]
    
    # The problem asks for the answer as a concatenated string of numbers.
    result_string = "".join(map(str, result))
    print(result_string)

solve_graphene_puzzle()