def solve():
    """
    Identifies and prints the incorrect statements about Attack Graphs.
    Based on the analysis:
    - Statement A is incorrect because LAG generation can be polynomial time, not necessarily exponential.
    - Statement D is incorrect because probabilities in cyclic graphs can be calculated, especially for SEGs which map to solvable Markov models.
    """
    incorrect_statements = ["A", "D"]
    
    # Sort the list alphabetically and join with a comma for the final output
    result = ",".join(sorted(incorrect_statements))
    
    print(result)

solve()