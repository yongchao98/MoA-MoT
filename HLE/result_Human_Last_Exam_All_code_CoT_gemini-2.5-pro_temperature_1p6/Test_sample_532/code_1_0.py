import math

def is_prime(n):
    """Helper function to check if a number is prime."""
    if n <= 1:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False
    for i in range(3, int(math.sqrt(n)) + 1, 2):
        if n % i == 0:
            return False
    return True

def check_if_filled_nilpotent(group_name: str):
    """
    Checks if a given group is a finite filled nilpotent group based on the known classification theorem.
    
    Args:
        group_name (str): A string representing the group.
                          Expected formats: "C_p" for cyclic group of prime order p (e.g., "C_5"),
                          "C3xC3" for C_3 x C_3, and "Q8" for the quaternion group.
    
    Prints the analysis and conclusion.
    """
    is_filled = False
    reason = f"'{group_name}' does not match the structure of a known finite filled nilpotent group."

    # Case 1: Cyclic group of prime order p, C_p
    if group_name.startswith("C_") and "_" in group_name:
        try:
            p_str = group_name.split("_")[1]
            if not p_str.isdigit():
                 raise ValueError("Non-integer value for p")
            p = int(p_str)
            
            # According to the theorem, C_n is filled iff n is a prime p
            # that satisfies the condition.
            if not is_prime(p):
                reason = f"The group C_{p} is not considered because {p} is not a prime number."
            else:
                # p is prime. Now we check the condition from the theorem.
                # Condition: p = 2 or p is a prime of the form 3k+2 (i.e., p % 3 == 2)
                if p == 2:
                    is_filled = True
                    reason = f"The group C_{p} is a filled nilpotent group because p=2."
                elif p % 3 == 2:
                    is_filled = True
                    reason = f"The group C_{p} is a filled nilpotent group because p={p} satisfies p % 3 == 2."
                else:
                    reason = f"The group C_{p} is not a filled nilpotent group because p={p} does not satisfy p=2 or p % 3 == 2."

        except (ValueError, IndexError):
            reason = f"Invalid format for cyclic group '{group_name}'. Use 'C_p' where p is a prime number."

    # Case 2: The group C_3 x C_3
    elif group_name == "C3xC3":
        is_filled = True
        reason = "The group C_3 x C_3 is a filled nilpotent group by classification."

    # Case 3: The quaternion group Q_8
    elif group_name == "Q8":
        is_filled = True
        reason = "The quaternion group Q_8 is a filled nilpotent group by classification."

    print(f"Group: {group_name} -> Is Filled Nilpotent: {is_filled}. Reason: {reason}")


if __name__ == '__main__':
    # A list of example groups to check against the classification.
    groups_to_test = [
        # Groups that are filled nilpotent
        "C_2",
        "C_5",     # 5 % 3 == 2
        "C_11",    # 11 % 3 == 2
        "C_17",    # 17 % 3 == 2
        "C3xC3",
        "Q8",
        
        # Groups that are NOT filled nilpotent
        "C_3",     # 3 % 3 == 0
        "C_4",     # 4 is not prime
        "C_6",     # 6 is not prime
        "C_7",     # 7 % 3 == 1
        "C_13",    # 13 % 3 == 1
        "D8",      # Dihedral group, not in the list
        "A5",      # Alternating group, not nilpotent
        "Invalid"
    ]
    
    print("--- Checking Groups Against the Classification of Finite Filled Nilpotent Groups ---")
    for group in groups_to_test:
        check_if_filled_nilpotent(group)