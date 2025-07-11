def solve_sequence():
    """
    This function solves a puzzle based on a known mathematical identity.
    The sequence given is one side of the smallest known solution to the Diophantine equation
    representing two sets of six integers whose fifth powers sum to the same total.
    The single integer that "completes" the sequence is this total sum.
    """
    
    # The sequence of numbers from the problem
    numbers = [24663, 35005, 119261, 196219, 211770, 227296]
    
    # Calculate the sum of the fifth powers of each number
    total_sum = sum(n**5 for n in numbers)
    
    # Build the equation string for the final output
    equation_parts = [f"{n}^5" for n in numbers]
    equation_str = " + ".join(equation_parts)
    
    print(f"{equation_str} = {total_sum}")

solve_sequence()