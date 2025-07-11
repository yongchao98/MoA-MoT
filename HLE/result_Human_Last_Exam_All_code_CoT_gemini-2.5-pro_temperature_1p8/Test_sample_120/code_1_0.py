import itertools

def solve_ita_matrix_code():
    """
    This function determines the correct ITA Matrix extension code based on a set of rules.
    It constructs the code for direct, Delta-operated flights eligible for
    Platinum Companion Certificates, formatted for minimal length and highest
    case-insensitive lexicographical value.
    """
    # 1. Define the components of the extension code.
    #    - o:DL specifies the operating carrier is Delta.
    #    - f bc=L|U|T|X|V filters for the required fare classes for the certificate.
    #    - f bc=0 filters for direct flights (0 connections).
    components = ['o:DL', 'f bc=L|U|T|X|V', 'f bc=0']

    # 2. Generate all possible orderings (permutations) of these components.
    permutations = list(itertools.permutations(components))

    # 3. Join the components in each permutation with a semicolon to form valid code strings.
    code_candidates = [";".join(p) for p in permutations]

    # 4. Sort the list of candidate strings case-insensitively.
    #    The last item in the sorted list will be the one with the highest lexicographical value.
    code_candidates.sort(key=lambda s: s.lower())
    
    # 5. The final answer is the last string in the sorted list.
    final_code = code_candidates[-1]
    
    print(final_code)

solve_ita_matrix_code()