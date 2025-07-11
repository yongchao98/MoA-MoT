import math

def solve():
    """
    Solves the problem for a given list of integers.
    The input is read from standard input.
    """
    import sys
    lines = sys.stdin.readlines()
    if not lines:
        # Handle empty input
        print("x = 1")
        print("total_length = 0")
        return

    a = [int(line.strip()) for line in lines if line.strip()]
    if not a:
        print("x = 1")
        print("total_length = 0")
        return
        
    a = [val for val in a if val > 0]
    if not a:
        print("x = 1")
        print("total_length = 0")
        return

    max_a = 0
    for val in a:
        if val > max_a:
            max_a = val

    # The delta array will store the difference S_q(x) - S_q(x-1)
    delta_S_q = [0] * (max_a + 2)

    for val in a:
        if val == 0:
            continue
        
        # Find all distinct quotients q = floor(val / k) > 0
        limit = int(math.sqrt(val))
        distinct_quotients = set()
        for k in range(1, limit + 1):
            distinct_quotients.add(val // k)
            distinct_quotients.add(k)
        
        # Sort quotients in descending order to process changes
        sorted_quotients = sorted(list(distinct_quotients), reverse=True)
        
        # The quotient for x > val is 0. The last positive quotient is 1.
        # The change from 1 to 0 happens at x = val + 1
        x_change = val + 1
        delta_S_q[x_change] += 0 - 1

        for q in sorted_quotients:
            if q == 0: continue
            # The quotient is q for x in (val/(q+1), val/q]
            # The change from a quotient > q to q happens at x = floor(val/q_prev) + 1
            # where q_prev is the quotient just larger than q.
            # A simpler way: The quotient drops from q at x = floor(val/q) + 1.
            x_change = (val // q) + 1
            if x_change > max_a + 1:
                continue

            new_q = val // x_change
            change = new_q - q
            if change != 0:
                delta_S_q[x_change] += change


    S_a = sum(a)
    S_q = S_a  # S_q(1) = sum(floor(a_i/1)) = sum(a_i)
    
    min_length = S_a # Length at x=1 is S_a + (1-1)*S_q = S_a
    best_x = 1

    for x in range(2, max_a + 2):
        S_q += delta_S_q[x]
        current_length = S_a + (1 - x) * S_q
        
        if current_length < min_length:
            min_length = current_length
            best_x = x

    print(f"The optimal value for x is: {best_x}")
    print(f"The minimized total length is: {min_length}")

solve()
