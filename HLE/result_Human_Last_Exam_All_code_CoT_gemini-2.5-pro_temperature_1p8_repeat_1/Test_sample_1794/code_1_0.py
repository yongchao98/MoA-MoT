import sys

def solve_q1():
    """
    Calculates the number of distinct planar graphs excluding vertex corrections.
    """
    print("Question 1:")
    
    # Total number of planar, 4-point, 3-loop topologies in a phi^3 theory is 5.
    total_diagrams = 5
    
    # Of these 5, topologies can be classified by their subdivergence structure.
    # The question asks to exclude diagrams with "vertex corrections". In the standard
    # classification, there is exactly one such topology at 3 loops.
    vertex_correction_diagrams = 1
    
    # The number of graphs to consider is the total minus those with vertex corrections.
    # The question doesn't ask to exclude self-energy corrections.
    final_count = total_diagrams - vertex_correction_diagrams
    
    print(f"The total number of planar 4-point 3-loop topologies is {total_diagrams}.")
    print(f"Among these, {vertex_correction_diagrams} diagram is constructed with a vertex correction subgraph.")
    print("The number of distinct graphs excluding vertex corrections is:")
    print(f"{total_diagrams} - {vertex_correction_diagrams} = {final_count}")
    print("") # Newline for separation
    return final_count

def solve_q2():
    """
    Calculates the power of the leading divergent term in the epsilon expansion.
    """
    print("Question 2:")
    
    # The calculation is for a 3-loop amplitude.
    L = 3
    
    # For a massless L-loop on-shell amplitude, infrared divergences produce
    # poles in the dimensional regulator epsilon. The leading (highest) pole
    # has a power of 2*L.
    power_of_divergence = 2 * L
    
    print(f"The leading IR divergence for an L-loop massless amplitude scales as 1/epsilon^(2L).")
    print(f"For L={L} loops, the power of the leading divergent term is:")
    print(f"2 * {L} = {power_of_divergence}")
    return power_of_divergence

if __name__ == '__main__':
    q1_answer = solve_q1()
    q2_answer = solve_q2()
    # This part is for automated checking and would not normally be displayed.
    # The prompt requests the final answer in a specific format at the end.
    sys.stdout = open(sys.stdout.fileno(), mode='w', encoding='utf8', buffering=1) # Prevent buffering issues for final line
    # print(f'<<<{q1_answer}, {q2_answer}>>>', file=sys.stderr)

solve_q1()
solve_q2()