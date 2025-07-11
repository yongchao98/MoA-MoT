def solve_arboricity_bounds():
    """
    This function determines the complexity classes for the arboricity of a subsampled graph.
    
    The analysis proceeds as follows:
    
    Case c = 1:
    - A lower bound is established using a graph of disjoint cliques. This shows that the arboricity can be as high as Omega(log(n)/log(log(n))).
    - This growth rate is omega(sqrt(log(n))) but o(log(n)), which corresponds to category 4.
    - An upper bound from literature shows the arboricity is O(log(n)).
    - Combining these, the tightest possible bound is not Theta(log(n)), making category 4 the best fit.
    - So, the first digit is 4.

    Case c = 2:
    - For c >= 2, the subsampled graph G' is known to have a very simple structure.
    - Specifically, with high probability, it does not contain any path of length 3.
    - This implies that every connected component is a star graph or a subgraph of one.
    - The arboricity of a star graph is 1.
    - Therefore, the arboricity of G' is O(1), which corresponds to category 1.
    - So, the second digit is 1.

    The final two-digit number is formed by concatenating the digits for f1 and f2.
    """
    
    # The first digit corresponds to the category for f1(n)
    f1_category = 4
    
    # The second digit corresponds to the category for f2(n)
    f2_category = 1
    
    # The final answer is the two-digit number
    final_answer = f1_category * 10 + f2_category
    
    # The problem asks to output each number in the final equation.
    # In this context, this means explaining the result.
    # Since the final output must be just the number, we print it directly.
    print(f"{f1_category}{f2_category}")

solve_arboricity_bounds()