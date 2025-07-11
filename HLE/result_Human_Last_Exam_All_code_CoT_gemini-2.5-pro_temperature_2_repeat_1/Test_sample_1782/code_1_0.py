def describe_partition(p):
    """Helper function to print a partition nicely."""
    return '{' + ', '.join(['{' + ', '.join(map(str, s)) + '}' for s in p]) + '}'

def main():
    """
    This function explores a finite analogue of the user's question.
    It demonstrates that for a finite sequence of refining partitions,
    a common refinement always exists.
    """
    N = 5
    base_set = set(range(N))
    print(f"Let's explore a finite version of the problem on the set S = {base_set}\n")

    # A sequence of partitions, where each is a refinement of the previous one.
    levels = []

    # Level 0: The whole set
    L0 = {frozenset(base_set)}
    levels.append(L0)
    print(f"Level 0 partition L_0 = {describe_partition(L0)}")

    # Level 1: Split the set into two
    L1 = {frozenset({0, 1}), frozenset({2, 3, 4})}
    levels.append(L1)
    print(f"Level 1 partition L_1 = {describe_partition(L1)} (refines L_0)")

    # Level 2: Refine the second part of L1
    L2 = {frozenset({0, 1}), frozenset({2}), frozenset({3, 4})}
    levels.append(L2)
    print(f"Level 2 partition L_2 = {describe_partition(L2)} (refines L_1)")

    # Level 3: Refine the first part of L2
    L3 = {frozenset({0}), frozenset({1}), frozenset({2}), frozenset({3, 4})}
    levels.append(L3)
    print(f"Level 3 partition L_3 = {describe_partition(L3)} (refines L_2)")

    # Level 4: The finest partition
    L4 = {frozenset({i}) for i in base_set}
    levels.append(L4)
    print(f"Level 4 partition L_4 = {describe_partition(L4)} (refines L_3)")

    print("\nWe have constructed a 'tree' of partitions of height 5.")

    # In a finite sequence of refinements, the last element refines all previous ones.
    common_refinement = levels[-1]

    print(f"\nDoes a common refinement of all levels exist?")
    print(f"Yes, one such common refinement is the last level itself: D = {describe_partition(common_refinement)}")

    print("\nVerification:")
    is_common_refinement = True
    for i, level in enumerate(levels):
        # Check if every set in D is a subset of some set in the current level.
        is_refinement = all(
            any(d.issubset(s) for s in level) for d in common_refinement
        )
        print(f"Is D a refinement of L_{i}? {is_refinement}")
        if not is_refinement:
            is_common_refinement = False

    print("\nConclusion for the finite case:")
    if is_common_refinement:
        print("A common refinement always exists for a finite sequence of refining partitions.")
        # Here we construct the "equation" showing the final result.
        # "Each number in the final equation" is interpreted as the elements of the sets.
        equation = f"Final_Answer_for_Finite_Case(S={describe_partition({base_set})}) = No"
        print("This means the answer to 'Does there always exist a tree without a common refinement?' is 'No' for finite sets.")
        print("\nFinal Result Illustrated:")
        print("L_4 refines L_3, which refines L_2, which refines L_1, which refines L_0.")
        # Print numbers from the "final equation" (i.e. the common refinement)
        final_elements = " ".join(sorted([item for subset in common_refinement for item in subset]))
        print(f"The numbers in the elements of the final partition D are: {final_elements}")


    else:
        # This case is not reachable in this setup
        print("A common refinement was not found.")

main()