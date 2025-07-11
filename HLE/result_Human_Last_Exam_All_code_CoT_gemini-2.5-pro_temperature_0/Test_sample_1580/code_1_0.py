import collections

def solve_logic_puzzle():
    """
    Solves the Aric and Pi number puzzle by simulating their reasoning process.
    """
    aric_num = 13
    pi_num = 10
    sums = [23, 20]

    # Step 1: Construct the chain of possible (Aric, Pi) pairs
    chain = collections.deque()
    
    # Build chain to the "left"
    curr_a, curr_p = aric_num, pi_num
    curr_sum_idx = 0 # Start with sum 23
    chain.append((curr_a, curr_p))
    
    # This loop generates the chain towards the (10,10) side
    a, p = curr_a, curr_p
    sum_idx = curr_sum_idx
    while True:
        # Pi's uncertainty
        other_sum = sums[1 - sum_idx]
        a_alt = other_sum - p
        if a_alt <= 0:
            break
        a = a_alt
        chain.append((a, p))
        sum_idx = 1 - sum_idx
        
        # Aric's uncertainty
        other_sum = sums[1 - sum_idx]
        p_alt = other_sum - a
        if p_alt <= 0:
            break
        p = p_alt
        chain.append((a, p))
        sum_idx = 1 - sum_idx

    # Build chain to the "right"
    # This loop generates the chain towards the (13,7) side
    a, p = curr_a, curr_p
    sum_idx = curr_sum_idx
    while True:
        # Aric's uncertainty
        other_sum = sums[1 - sum_idx]
        p_alt = other_sum - a
        if p_alt <= 0:
            break
        p = p_alt
        chain.appendleft((a, p))
        sum_idx = 1 - sum_idx

        # Pi's uncertainty
        other_sum = sums[1 - sum_idx]
        a_alt = other_sum - p
        if a_alt <= 0:
            break
        a = a_alt
        chain.appendleft((a, p))
        sum_idx = 1 - sum_idx

    # Step 2: Simulate the game
    day = 1
    while True:
        # Aric's turn
        # Aric can solve if his number is at the left endpoint of the chain
        # because his alternative has been eliminated by the previous player's pass.
        # The very first turn is a special case where the alternative is impossible (non-positive).
        left_endpoint_a, left_endpoint_p = chain[0]
        if aric_num == left_endpoint_a:
            winner = "Aric"
            reason_p = sums[0] - aric_num if (aric_num + left_endpoint_p) == sums[0] else sums[1] - aric_num
            final_sum = aric_num + pi_num
            break
        
        # Aric passes, eliminating the left endpoint
        chain.popleft()

        # Pi's turn
        # Pi can solve if her number is at the right endpoint
        right_endpoint_a, right_endpoint_p = chain[-1]
        if pi_num == right_endpoint_p:
            winner = "Pi"
            reason_a = sums[0] - pi_num if (right_endpoint_a + pi_num) == sums[0] else sums[1] - pi_num
            final_sum = aric_num + pi_num
            break

        # Pi passes, eliminating the right endpoint
        chain.pop()
        
        day += 1

    # Step 3: Print the results
    print(f"{winner} gives the answer on Day {day}.")
    print(f"Here is {winner}'s reasoning:")
    if winner == "Aric":
        print(f"\"My number is {aric_num}. Pi's number could be {pi_num} (sum {sums[0]}) or {reason_p} (sum {sums[1]}).\"")
        print(f"\"However, I know that on Day {day-1}, Pi passed. If her number had been {reason_p}, she would have been able to determine the sum because her other possibility was eliminated by my pass on that same day.\"")
        print(f"\"Since she passed, her number cannot be {reason_p}. Therefore, her number must be {pi_num}.\"")
    else: # Pi's reasoning
        print(f"\"My number is {pi_num}. Aric's number could be {aric_num} (sum {sums[0]}) or {reason_a} (sum {sums[1]}).\"")
        print(f"\"However, I know that on Day {day}, Aric passed. If his number had been {reason_a}, he would have been able to determine the sum because his other possibility was eliminated by my pass on Day {day-1}.\"")
        print(f"\"Since he passed, his number cannot be {reason_a}. Therefore, his number must be {aric_num}.\"")
    
    print("\nThe final equation is:")
    print(f"{aric_num} + {pi_num} = {final_sum}")


solve_logic_puzzle()