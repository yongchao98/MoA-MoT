def check_answer():
    """
    Checks the correctness of the answer to the chemistry question.
    A compound is optically active if it is chiral and not a meso compound.
    """

    # The question's options are A) 2, B) 4, C) 5, D) 3.
    # The given answer is 'D', which corresponds to a count of 3.
    expected_count = 3
    
    # Analysis of each compound based on chemical principles of chirality.
    # is_active: True if the compound is chiral and optically active.
    # is_active: False if the compound is achiral (e.g., has a plane of symmetry, is meso).
    compounds_analysis = [
        {
            "name": "(Z)-1-chloro-2-methylbut-1-ene",
            "is_active": False,
            "reason": "Achiral. Has a plane of symmetry (planar around the double bond)."
        },
        {
            "name": "(3aR,7aS,E)-8-(chloromethylene)hexahydro-4,7-methanoisobenzofuran-1,3-dione",
            "is_active": True,
            "reason": "Chiral. The name specifies a single enantiomer of a complex, asymmetric molecule."
        },
        {
            "name": "(2R,3S)-2,3-dimethylsuccinic acid",
            "is_active": False,
            "reason": "Achiral. It is a meso compound with an internal plane of symmetry."
        },
        {
            "name": "(2R,3R)-2,3-dimethylsuccinic acid",
            "is_active": True,
            "reason": "Chiral. This diastereomer lacks a plane of symmetry."
        },
        {
            "name": "(R)-cyclohex-3-en-1-ol",
            "is_active": True,
            "reason": "Chiral. The name specifies a single enantiomer with a chiral center at C1."
        },
        {
            "name": "(1s,3s,5s)-cyclohexane-1,3,5-triol",
            "is_active": False,
            "reason": "Achiral. The all-cis configuration results in a highly symmetric molecule with planes of symmetry."
        },
        {
            "name": "1-cyclopentyl-3-methylbutan-1-one",
            "is_active": False,
            "reason": "Achiral. No chiral centers and has a plane of symmetry."
        }
    ]

    # Calculate the number of optically active compounds from the analysis.
    calculated_count = sum(1 for compound in compounds_analysis if compound["is_active"])

    # Check if the calculated count matches the expected count from the answer.
    if calculated_count == expected_count:
        return "Correct"
    else:
        active_compounds = [c["name"] for c in compounds_analysis if c["is_active"]]
        inactive_compounds = [c["name"] for c in compounds_analysis if not c["is_active"]]
        
        reason = f"The answer is incorrect. "
        reason += f"The provided answer corresponds to {expected_count} optically active compounds, but the correct analysis shows there are {calculated_count}. "
        reason += f"The optically active compounds are: {active_compounds}. "
        reason += f"The optically inactive compounds are: {inactive_compounds}."
        return reason

# Run the check
print(check_answer())