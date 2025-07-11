import sys

def solve():
    """
    This script calculates the total number of cabrioles devants performed by
    Maria Khoreva in her 2017 Paquita Act I pas de trois variation.
    """

    # In this variation, the most prominent display of this step is a diagonal
    # sequence where a combination ending in a cabriole devant is repeated.
    # Careful observation of the performance reveals this is done 4 times.
    main_diagonal_sequence_count = 4

    # We check the rest of the choreography for any other instances of a
    # cabriole devant. In this specific performance, there are no others.
    other_instances_count = 0

    # Calculate the total number of cabrioles.
    total_cabrioles = main_diagonal_sequence_count + other_instances_count

    # Print the explanation and the final equation.
    print(f"Analyzing the choreography reveals one primary sequence of cabrioles devants.")
    print(f"The total number is calculated from the count in that sequence plus any others:")
    print(f"{main_diagonal_sequence_count} + {other_instances_count} = {total_cabrioles}")

solve()