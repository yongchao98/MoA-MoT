import math

def calculate_bottcher_complexity():
    """
    Calculates the Böttcher Molecular Complexity for the product of the
    Favorskii rearrangement of 2-chlorocyclohexanone.
    """
    # Step 1: Identify the product and represent its molecular graph.
    # The product is cyclopentanecarboxylic acid. Its non-hydrogen skeleton is
    # represented here by an adjacency list, where keys are atoms and values
    # are the atoms they are bonded to. Multiple bonds are treated as single edges.
    molecular_graph = {
        'C_ring_1': ['C_ring_2', 'C_ring_5', 'C_carboxyl'],
        'C_ring_2': ['C_ring_1', 'C_ring_3'],
        'C_ring_3': ['C_ring_2', 'C_ring_4'],
        'C_ring_4': ['C_ring_3', 'C_ring_5'],
        'C_ring_5': ['C_ring_4', 'C_ring_1'],
        'C_carboxyl': ['C_ring_1', 'O_carbonyl', 'O_hydroxyl'],
        'O_carbonyl': ['C_carboxyl'],
        'O_hydroxyl': ['C_carboxyl']
    }
    
    # Step 2: Determine the parameters for the calculation.
    
    # Calculate degrees (delta_i) for each non-hydrogen atom.
    # The degree is the number of connections for each atom in the graph.
    degrees = {atom: len(neighbors) for atom, neighbors in molecular_graph.items()}
    
    # Calculate N_b, the total number of bonds (edges) in the graph.
    # This is the sum of all degrees divided by 2.
    N_b = sum(degrees.values()) / 2
    N_b = int(N_b)
    
    # Calculate the sum of the squares of the degrees.
    sum_of_squares_degrees = sum(d**2 for d in degrees.values())
    
    # Step 3: Apply the Böttcher formula and print the results.
    print("The product of the Favorskii rearrangement of 2-chlorocyclohexanone is cyclopentanecarboxylic acid.")
    print("\nCalculating the Böttcher Molecular Complexity:")
    print("-" * 45)
    
    print(f"1. Number of bonds (N_b) between non-hydrogen atoms: {N_b}")
    
    sorted_degrees = sorted(degrees.values(), reverse=True)
    print(f"2. Degrees of the 8 non-hydrogen atoms (delta_i): {sorted_degrees}")
    
    print(f"3. Sum of squares of degrees (sum(delta_i^2)): {sum_of_squares_degrees}")
    
    N_b_squared = N_b**2
    numerator = N_b_squared - sum_of_squares_degrees
    
    print("\nUsing the formula: BC = (N_b^2 - sum(delta_i^2)) / N_b^2")
    print(f"BC = ({N_b}^2 - {sum_of_squares_degrees}) / {N_b}^2")
    print(f"BC = ({N_b_squared} - {sum_of_squares_degrees}) / {N_b_squared}")
    print(f"BC = {numerator} / {N_b_squared}")
    
    # Final calculation
    bottcher_complexity = numerator / N_b_squared
    
    print(f"\nThe Böttcher Molecular Complexity is: {bottcher_complexity}")
    
    return bottcher_complexity

# Execute the function
final_answer = calculate_bottcher_complexity()
# The required output format is now outside the function.
# This makes it clearer how the final answer is produced.
# print(f"\n<<<{final_answer}>>>")

if __name__ == '__main__':
    pass
