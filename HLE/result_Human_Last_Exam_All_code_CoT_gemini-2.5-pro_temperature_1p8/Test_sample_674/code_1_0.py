def solve_sequence():
    """
    Deduces and computes the sequence S4.
    The rule is deduced to be a Hofstadter-Conway-like sequence:
    s[n] = s[s[n-1]] + s[n-s[n-1]]
    """
    s = {1: 1, 2: 1}
    num_terms = 42

    print(f"s[1] = 1")
    print(f"s[2] = 1")

    for n in range(3, num_terms + 1):
        try:
            prev_s1 = s[n-1]
            
            # The indices for the recurrence relation
            idx1 = prev_s1
            idx2 = n - prev_s1

            if idx1 not in s or idx2 not in s:
                print(f"Cannot compute s[{n}]: missing required terms s[{idx1}] or s[{idx2}]")
                break

            val1 = s[idx1]
            val2 = s[idx2]
            
            result = val1 + val2
            s[n] = result
            
            # Print the equation for the current step
            print(f"s[{n}] = s[s[{n-1}]] + s[{n}-s[{n-1}]] = s[{idx1}] + s[{idx2}] = {val1} + {val2} = {result}")

        except KeyError:
            print(f"Error: Required term for s[{n}] not found in the sequence.")
            break
            
solve_sequence()
# The final answer is the formula itself.
print("\n<<<R(s[n]) = s[s[n-1]] + s[n - s[n-1]]>>>")
