import itertools

def solve_isomer_problem():
    """
    Calculates the number of isomers for the dinuclear complex
    [(bpy)2Ru-(dptztz)-Ru(bpy)2]^4+.
    """
    print("Step 1: Identifying the building blocks for the dinuclear complex.")
    print("Each Ru center is chiral (Δ or Λ) and binds an unsymmetrical chelate from the bridge,")
    print("leading to two diastereomers for each chirality. This gives 4 building blocks:\n")

    blocks = ["ΔA", "ΔB", "ΛA", "ΛB"]
    enantiomers = {"ΔA": "ΛA", "ΛA": "ΔA", "ΔB": "ΛB", "ΛB": "ΔB"}

    print(f"Building Blocks: {blocks}")
    print("Enantiomeric pairs are (ΔA, ΛA) and (ΔB, ΛB).\n")


    print("Step 2: Combining building blocks to form dinuclear complexes of the type {Ru1, Ru2}.")
    print("Since the bridging ligand is symmetric, {Ru1, Ru2} is the same as {Ru2, Ru1}.\n")

    # Generate all unique pairs (combinations with replacement)
    all_pairs = list(itertools.combinations_with_replacement(blocks, 2))

    print(f"All possible unique combinations: {all_pairs}\n")

    print("Step 3: Classifying each combination as meso or part of a chiral pair and counting.\n")

    found_isomers = set()
    total_isomers = 0
    sum_equation = []

    for pair in all_pairs:
        # Sort to handle commutativity, e.g. ('ΔA', 'ΛA') is same as ('ΛA', 'ΔA')
        p1, p2 = sorted(pair)

        # If this pair or its enantiomer has already been processed, skip.
        if (p1, p2) in found_isomers:
            continue

        enantiomer_p1 = enantiomers[p1]
        enantiomer_p2 = enantiomers[p2]
        enantiomer_pair = tuple(sorted((enantiomer_p1, enantiomer_p2)))

        # Check for meso compound
        if (p1, p2) == enantiomer_pair:
            print(f"Combination {{{p1}, {p2}}} is a meso compound.")
            print(f"  - Mirror image is {{{enantiomer_p1}, {enantiomer_p2}}}, which is the same molecule.")
            print("  - This counts as 1 isomer.")
            found_isomers.add((p1, p2))
            total_isomers += 1
            sum_equation.append("1")
        else:
            # It's a chiral pair
            print(f"Combination {{{p1}, {p2}}} is part of a chiral pair.")
            print(f"  - Its enantiomer is {{{enantiomer_pair[0]}, {enantiomer_pair[1]}}}.")
            print("  - This counts as 2 isomers.")
            found_isomers.add((p1, p2))
            found_isomers.add(enantiomer_pair)
            total_isomers += 2
            sum_equation.append("2")
        print("-" * 20)

    print("\nFinal Calculation:")
    final_equation = " + ".join(sum_equation)
    print(f"The total number of isomers is the sum of isomers from each unique combination:")
    print(f"{final_equation} = {total_isomers}")

solve_isomer_problem()
<<<10>>>