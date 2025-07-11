import collections

def solve_geology_quiz():
    """
    Analyzes geological statements and reports the consensus/debated status.
    The phrase 'output each number in the final equation' from the prompt is interpreted
    as showing the classification for each numbered statement.
    """

    # A dictionary mapping statement number to its classification and justification.
    # 'C' for Consensus, 'D' for Debated.
    statements = {
        1: ("C", "The Morrison formation as a foredeep/foreland basin deposit is a consensus view."),
        2: ("D", "The cause of metamorphic core complexes is highly debated, with the slab window model being one of several competing hypotheses."),
        3: ("D", "The 'Nevadaplano' is an influential but debated model; the elevation and extent of the proposed plateau are not settled."),
        4: ("D", "The 'Arizonaplano' and its proposed formation via magmatic thickening is a specific, debated hypothesis."),
        5: ("C", "The spatial relationship (inboard Laramide, outboard Sevier) is a fundamental, consensus observation."),
        6: ("C", "The origin of the Sierra Nevada and Idaho batholiths from the ancestral Cordilleran arc is a foundational concept."),
        7: ("C", "The southwestward propagation of volcanism (the Yellowstone hotspot track) is a well-documented, consensus model."),
        8: ("D", "The uplift history of the Colorado Plateau is a classic, long-standing geological debate."),
        9: ("D", "The detailed uplift history of the Cascades range is complex and debated."),
        10: ("D", "The Shatsky conjugate as the specific driver for the Laramide orogeny is a debated hypothesis within the broader flat-slab subduction model.")
    }
    
    # Sort statements by number to ensure correct order
    sorted_statements = collections.OrderedDict(sorted(statements.items()))
    
    final_answer_list = []
    
    print("Classification for each statement:")
    for num, (classification, reason) in sorted_statements.items():
        print(f"({num}) -> {classification}")
        final_answer_list.append(classification)
        
    final_string = "".join(final_answer_list)
    print("\nFinal combined string:")
    print(final_string)

solve_geology_quiz()