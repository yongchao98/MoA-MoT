import itertools

def find_optimal_extension_code():
    """
    This script determines the optimal ITA Matrix extension code based on a set of rules.

    Rules:
    1.  Direct (non-stop) Delta flights.
    2.  Eligible for Platinum-tier Delta Companion Certificates (Fare classes L, T, U, V, X).
    3.  Minimize overall string length.
    4.  Among strings of minimal length, choose the one with the highest case-insensitive lexicographical value.
    """

    # Components of the extension code
    # 1. Airline: Delta
    airline_code = "DL"

    # 2. Direct Flight (non-stop): f ~1 is the most compact syntax for one segment.
    non_stop_code = "f ~1"

    # 3. Fare Classes: Platinum Companion Certificate classes (L,U,T,X,V), alphabetized and lowercased.
    # The 'v' flag specifies booking class.
    fare_class_code = "f v l|t|u|v|x"

    # List of components to be arranged
    components = [fare_class_code, non_stop_code, airline_code]

    # Generate all possible orderings (permutations) of the components
    permutations = list(itertools.permutations(components))

    # Join the components with a single space to create the full code strings
    possible_codes = [" ".join(p) for p in permutations]

    # Sort the list of codes lexicographically (case-insensitively) in ascending order
    # The last element in the sorted list will be the one with the highest value.
    # str.casefold() is used for robust case-insensitive comparison.
    sorted_codes = sorted(possible_codes, key=str.casefold)

    # The optimal code is the last one in the sorted list
    optimal_code = sorted_codes[-1]

    print("The required ITA Matrix outbound extension code is:")
    print(optimal_code)

find_optimal_extension_code()