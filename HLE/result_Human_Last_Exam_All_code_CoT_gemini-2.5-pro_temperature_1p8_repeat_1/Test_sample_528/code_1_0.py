import math

def solve():
    """
    This function calculates the number of sets (S1, S2, S3, S4) satisfying
    S1 ⊂ S2 ⊂ S3 ⊂ S4 ⊂ {1,2,3,4,5} and i ∈ Si for i = 1,2,3.
    """

    # The universal set is U = {1, 2, 3, 4, 5}.
    # The chain of subsets S1 ⊂ S2 ⊂ S3 ⊂ S4 ⊂ U partitions U into 5 disjoint regions.
    # We count the number of possible regions for each element.

    # Choices for element 1: 1 ∈ S1. It must be in the region S1.
    choices_1 = 1

    # Choices for element 2: 2 ∈ S2. It can be in S1 or S2 \ S1.
    choices_2 = 2

    # Choices for element 3: 3 ∈ S3. It can be in S1, S2 \ S1, or S3 \ S2.
    choices_3 = 3

    # Choices for element 4: No restrictions. It can be in any of the 5 regions.
    # (S1, S2 \ S1, S3 \ S2, S4 \ S3, U \ S4)
    choices_4 = 5

    # Choices for element 5: No restrictions. Same as element 4.
    choices_5 = 5

    # The total number of ways is the product of the choices for each element.
    total_sets = choices_1 * choices_2 * choices_3 * choices_4 * choices_5

    print(f"The number of sets is the product of the number of choices for each element based on the given constraints.")
    print(f"Choices for element 1 (must be in S1): {choices_1}")
    print(f"Choices for element 2 (must be in S2): {choices_2}")
    print(f"Choices for element 3 (must be in S3): {choices_3}")
    print(f"Choices for element 4 (no restrictions): {choices_4}")
    print(f"Choices for element 5 (no restrictions): {choices_5}")
    print("\nFinal equation:")
    print(f"{choices_1} * {choices_2} * {choices_3} * {choices_4} * {choices_5} = {total_sets}")
    print(f"\nThus, there are {total_sets} such sets.")


solve()