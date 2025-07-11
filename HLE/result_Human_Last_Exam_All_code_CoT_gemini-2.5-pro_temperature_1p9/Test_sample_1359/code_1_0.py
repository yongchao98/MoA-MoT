def solve_missing_triplets():
    """
    Calculates the missing triplets and their sum based on the provided rules.

    The solution proceeds based on the assumption that the Horizontal Transformation
    rules should be used to complete the third row, as the provided data is
    inconsistent with a combined horizontal and vertical rule application.
    """

    # Starting triplet for Row 3
    t1 = [7, 2, 9]
    print(f"The first known triplet in the third row is {t1}.\n")

    missing_elements = []

    # --- Calculate the second triplet (T2) from the first (T1) ---
    print("--- Calculating the middle triplet [? ? ?] from [7 2 9] ---")
    x1, y1, z1 = t1
    
    # Check the horizontal transformation condition
    condition_sum = x1 + y1
    print(f"Checking condition: x + y = {x1} + {y1} = {condition_sum}")

    if condition_sum > 10:
        # This block will not be executed for T1 -> T2
        print("Condition x + y > 10 met.")
        x2 = (x1 * 3 - y1) % 12
        y2 = (y1 * 2 + 4) % 12
        z2 = (z1 + x1) % 12
    else:
        print("Condition x + y <= 10 met.")
        print(f"Next x = ({x1} * 2 + {y1}) mod 12")
        x2 = (x1 * 2 + y1) % 12
        print(f"         = ({x1*2} + {y1}) mod 12 = {x1*2+y1} mod 12 = {x2}")

        print(f"Next y = ({y1} * 3 - 2) mod 12")
        y2 = (y1 * 3 - 2) % 12
        print(f"         = ({y1*3} - 2) mod 12 = {y1*3-2} mod 12 = {y2}")

        print(f"Next z = ({z1} * 2) mod 12")
        z2 = (z1 * 2) % 12
        print(f"         = {z1*2} mod 12 = {z2}")

    t2 = [x2, y2, z2]
    missing_elements.extend(t2)
    print(f"\nThe calculated middle triplet is {t2}.\n")

    # --- Calculate the third triplet (T3) from the second (T2) ---
    print(f"--- Calculating the right triplet [? ? ?] from {t2} ---")
    
    # Check the horizontal transformation condition for T2
    condition_sum = x2 + y2
    print(f"Checking condition: x + y = {x2} + {y2} = {condition_sum}")

    if condition_sum > 10:
        # This block will not be executed for T2 -> T3
        print("Condition x + y > 10 met.")
        x3 = (x2 * 3 - y2) % 12
        y3 = (y2 * 2 + 4) % 12
        z3 = (z2 + x2) % 12
    else:
        print("Condition x + y <= 10 met.")
        print(f"Next x = ({x2} * 2 + {y2}) mod 12")
        x3 = (x2 * 2 + y2) % 12
        print(f"         = ({x2*2} + {y2}) mod 12 = {x2*2+y2} mod 12 = {x3}")
        
        print(f"Next y = ({y2} * 3 - 2) mod 12")
        y3 = (y2 * 3 - 2) % 12
        print(f"         = ({y2*3} - 2) mod 12 = {y2*3-2} mod 12 = {y3}")
        
        print(f"Next z = ({z2} * 2) mod 12")
        z3 = (z2 * 2) % 12
        print(f"         = {z2*2} mod 12 = {z3}")
    
    t3 = [x3, y3, z3]
    missing_elements.extend(t3)
    print(f"\nThe calculated right triplet is {t3}.\n")

    # --- Sum the missing elements ---
    total_sum = sum(missing_elements)
    print("--- Calculating the final sum ---")
    print(f"The missing elements are: {missing_elements[0]}, {missing_elements[1]}, {missing_elements[2]}, {missing_elements[3]}, {missing_elements[4]}, {missing_elements[5]}")
    sum_string = " + ".join(map(str, missing_elements))
    print(f"Sum = {sum_string} = {total_sum}")

solve_missing_triplets()