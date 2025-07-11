import sys

def solve_geopolitical_graph_problem():
    """
    Solves the problem regarding the chromatic number and planarity of the Asian country graph.
    """
    # Step 1: Determine the chromatic number of the Asian graph BEFORE the Soviet dissolution.
    # A map of countries is a planar graph. By the Four-Color Theorem, its chromatic number is at most 4.
    # Before 1991, the map of Asia, simplified by the large and contiguous USSR, did not possess
    # an obvious subgraph (like a 4-clique or certain wheel graphs) that would necessitate 4 colors.
    # The lack of such a forcing structure suggests the map was 3-colorable.
    chromatic_number_before = 3
    print(f"Chromatic number of the Asian subgraph before Soviet dissolution: {chromatic_number_before}")

    # Step 2: Determine the chromatic number of the Asian graph AFTER the Soviet dissolution.
    # The dissolution created new countries, introducing a subgraph that requires 4 colors.
    # This structure is centered on Uzbekistan and its neighbors: Kazakhstan, Kyrgyzstan,
    # Tajikistan, Afghanistan, and Turkmenistan. These 5 neighbors form a 5-sided cycle (an odd cycle)
    # on the "rim" of the wheel, which requires 3 colors. Uzbekistan, the "hub", is adjacent to all
    # 5 countries on the rim, forcing the use of a 4th color.
    chromatic_number_after = 4
    print(f"Chromatic number of the Asian subgraph after Soviet dissolution: {chromatic_number_after}")

    # Step 3: Calculate delta_soviet, the incremental change.
    delta_soviet = chromatic_number_after - chromatic_number_before
    print(f"The incremental change in the chromatic number, delta_soviet, is {chromatic_number_after} - {chromatic_number_before} = {delta_soviet}")

    # Step 4: Determine beta, representing the change in planarity.
    # A graph of contiguous country borders is inherently planar.
    # The graph was planar before 1991 and remained planar after. The planarity did not change.
    # As per the problem, if planarity does not change, beta is 1.
    beta = 1
    print(f"The planarity of the graph did not change, therefore beta = {beta}")

    # Step 5: Calculate the final answer.
    final_answer = beta * delta_soviet
    print(f"The final answer is beta * delta_soviet = {beta} * {delta_soviet} = {final_answer}")
    
    # Required final output format
    sys.stdout.write(f"<<<{final_answer}>>>")

solve_geopolitical_graph_problem()