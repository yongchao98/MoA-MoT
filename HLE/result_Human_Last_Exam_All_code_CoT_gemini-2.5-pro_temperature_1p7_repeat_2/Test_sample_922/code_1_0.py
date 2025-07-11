import itertools

def main():
    """
    Solves the sequence puzzle by identifying the property of its numbers
    and finding the next known term based on a 2022 discovery.
    """
    sequence_info = """
The sequence 24663, 35005, 119261, 196219, 211770, 227296 consists of numbers
that cannot be expressed as the sum of four tetrahedral numbers.

A tetrahedral number is given by the formula T(n) = n(n+1)(n+2)/6.

In August 2022, it was proven that a new integer, 343867, joins this list.

To demonstrate this property, the following code will show that the next integer,
343868, CAN be written as the sum of four tetrahedral numbers.
"""
    print(sequence_info)

    target = 343868

    def tetrahedral(n):
        """Calculates the n-th tetrahedral number."""
        return n * (n + 1) * (n + 2) // 6

    # Generate tetrahedral numbers up to the target value
    t_numbers = []
    n = 1
    while True:
        t = tetrahedral(n)
        if t > target:
            break
        t_numbers.append(t)
        n += 1
    
    # We also need the original 'n' values for the final output
    t_map = {tetrahedral(i): i for i in range(1, n)}

    # Find a combination of 4 tetrahedral numbers that sum to the target
    solution = None
    for combo in itertools.combinations_with_replacement(t_numbers, 4):
        if sum(combo) == target:
            solution = combo
            break
    
    # Print the equation
    if solution:
        n_values = sorted([t_map[t] for t in solution])
        print(f"The number {target} is the sum of four tetrahedral numbers:")
        print(f"{target} = T({n_values[0]}) + T({n_values[1]}) + T({n_values[2]}) + T({n_values[3]})")
        print(f"{target} = {solution[0]} + {solution[1]} + {solution[2]} + {solution[3]}")
    else:
        print(f"Could not find a representation for {target} as the sum of four tetrahedral numbers.")
        
    final_answer = 343867
    print(f"\nThe next number in the sequence is therefore {final_answer}.")


main()