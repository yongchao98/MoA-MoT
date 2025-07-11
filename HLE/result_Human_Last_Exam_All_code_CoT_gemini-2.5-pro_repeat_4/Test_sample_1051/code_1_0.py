import sys
import math
import bisect

def solve():
    """
    Solves the optimization problem to find the integer x that minimizes the total length.
    """
    try:
        lines = sys.stdin.readlines()
        if not lines:
            # Handle empty input
            print("best x: 1")
            print("minimized length: 0")
            return
            
        a_list = [int(line.strip()) for line in lines if line.strip()]
        if not a_list:
            # Handle input with only empty lines
            print("best x: 1")
            print("minimized length: 0")
            return

        n = len(a_list)
        total_sum_a = sum(a_list)
        max_a = 0
        if a_list:
            max_a = max(a_list)

        # 1. Generate candidate points
        candidates = {1}
        if max_a > 0:
            for a in a_list:
                if a == 0: continue
                limit = int(math.sqrt(a))
                for k in range(1, limit + 1):
                    candidates.add(k)
                    candidates.add(a // k)
        
        # We also need to check x-1 for each candidate x, as the length function can dip between integers
        # Adding x-1 for each candidate where x>1 is a safer way to probe the function's behavior
        # However, the logic that the minimum of the piecewise linear function L(x) must be at the
        # endpoints of the constant-F(x) intervals is sound. Those endpoints are the candidates.
        # So we just test the candidates themselves.

        sorted_candidates = sorted(list(candidates))
        m = len(sorted_candidates)
        
        # 2. Efficiently compute F(p) for all p in sorted_candidates
        # F(p) = sum(floor(a_i / p))
        # We use a difference array for efficient range updates.
        diff_array = [0] * (m + 1)

        for a in a_list:
            if a == 0: continue
            s = int(math.sqrt(a))
            
            # Partition the domain of x into [1, s] and (s, a]
            # Processing x in [1, s]
            last_k = a # k for x=1 is a
            # For x = 1, quotient is 'a'
            idx_1 = bisect.bisect_left(sorted_candidates, 1)
            diff_array[idx_1] += a
            if idx_1 + 1 <= m:
                diff_array[idx_1 + 1] -= a

            # For x in [2, s]
            for x in range(2, s + 1):
                k = a // x
                # Add k for a single point x
                idx = bisect.bisect_left(sorted_candidates, x)
                if idx < m and sorted_candidates[idx] == x:
                    diff_array[idx] += k
                    if idx + 1 <= m:
                        diff_array[idx+1] -= k
            
            # Processing x > s. This corresponds to quotients k < sqrt(a)
            for k in range(1, s + 1):
                x_max = a // k
                x_min = max(s + 1, a // (k + 1) + 1)
                
                if x_min > x_max:
                    continue

                idx_start = bisect.bisect_left(sorted_candidates, x_min)
                idx_end = bisect.bisect_right(sorted_candidates, x_max) - 1

                if idx_start <= idx_end:
                    diff_array[idx_start] += k
                    if idx_end + 1 <= m:
                        diff_array[idx_end + 1] -= k

        # 3. Reconstruct F_values from the difference array
        f_values = [0] * m
        current_f = 0
        for i in range(m):
            current_f += diff_array[i]
            f_values[i] = current_f

        # 4. Calculate L(p) for all p and find the minimum
        min_len = float('inf')
        best_x = 1

        for i in range(m):
            p = sorted_candidates[i]
            if p == 0: continue # x must be positive
            
            # L(p) = sum(a_i) + (1-p) * F(p)
            current_len = total_sum_a + (1 - p) * f_values[i]

            if current_len < min_len:
                min_len = current_len
                best_x = p
        
        # The problem is defined for positive integers x, so min_len should be an integer
        min_len = round(min_len)

        print(f"The positive integer x that minimizes the total length is: {best_x}")
        print(f"The minimized total length is: {min_len}")

    except Exception as e:
        # On error, print a message.
        print(f"An error occurred: {e}", file=sys.stderr)

solve()
