import itertools

def solve_ita_matrix_code():
    """
    This script determines the correct ITA Matrix extension code based on a specific set of rules.
    The rules are:
    1. Airline: Delta (DL)
    2. Stops: Direct (no layovers)
    3. Fare Classes: Eligible for Platinum Delta Companion Certificate (L, U, T, X, V)
    4. Optimization: Minimal string length, then highest case-insensitive lexicographical sort.
    """

    # Step 1: Define the components of the extension code in their shortest form.
    # Airline component
    airline_code = "DL"

    # Direct flight component (shortest form)
    stops_code = "stops 0"

    # Fare class component for Platinum Companion Certificate (L,U,T,X,V)
    # The format 'f bc=A|B|C' is the most compact.
    fare_class_code = "f bc=L|U|T|X|V"

    components = [airline_code, stops_code, fare_class_code]

    # Step 2: Generate all possible orderings (permutations) of the components.
    all_permutations = list(itertools.permutations(components))

    # Step 3: Format each permutation into a full extension code string.
    # The standard separator is a semicolon followed by a space.
    candidate_codes = [" ; ".join(p) for p in all_permutations]

    # Step 4: Find the code that is lexicographically highest.
    # The max() function on a list of strings performs this comparison.
    final_code = max(candidate_codes)

    # Step 5: Print the final result.
    print("The ITA Matrix extension code is:")
    print(final_code)

solve_ita_matrix_code()