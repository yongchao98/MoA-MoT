def check_cardinal_property(kappa_name: str):
    """
    Checks if a cardinal kappa is regular or singular and explains the implication
    for the set theory problem.
    The input string should be of the form 'aleph_n' where n is an integer or 'omega'.
    """
    
    print(f"--- Checking for kappa = {kappa_name} ---")
    
    is_regular = False
    reason = ""
    
    if kappa_name.startswith("aleph_"):
        index_str = kappa_name.split("_")[1]
        
        if index_str.isdigit():
            # For any integer n, aleph_n is a regular cardinal.
            # aleph_0 is regular by definition.
            # For n > 0, n is a successor ordinal, and aleph_n is a successor cardinal,
            # which is always regular.
            is_regular = True
            reason = f"{kappa_name} is a regular cardinal."
        elif index_str == "omega":
            # aleph_omega is the first singular cardinal.
            # It's the supremum of a countable sequence of smaller cardinals.
            # cf(aleph_omega) = omega < aleph_omega.
            is_regular = False
            reason = f"{kappa_name} is a singular cardinal."
        else:
            print(f"Cannot determine regularity for {kappa_name}.")
            return
    else:
        print(f"Unrecognized cardinal format: {kappa_name}.")
        return

    print(reason)
    
    if is_regular:
        print("Result: A function f with the specified property EXISTS for this kappa.")
    else:
        print("Result: A function f with the specified property does NOT exist for this kappa.")
    print("")

def main():
    # Test with aleph_0 = omega, which is a regular cardinal.
    check_cardinal_property("aleph_0")
    
    # Test with aleph_5, which is a successor cardinal and thus regular.
    check_cardinal_property("aleph_5")

    # Test with aleph_omega, which is the first singular cardinal.
    check_cardinal_property("aleph_omega")

if __name__ == "__main__":
    main()