def solve():
    """
    Finds the tuple (a, b, c, d) that maximizes the Ducci sequence length
    and has the smallest sum, then computes (a + b - c - d) mod 1000.
    """
    limit = 10_000_000

    # Generate the Tribonacci sequence starting with (0, 0, 1) as it grows the slowest,
    # allowing for a higher n (and thus longer f) within the limit.
    trib = [0, 0, 1]
    while trib[-1] <= limit:
        next_trib = trib[-1] + trib[-2] + trib[-3]
        trib.append(next_trib)

    # The last element is over the limit, so we don't use it.
    # The tuple (a, b, c, d) that maximizes f is formed by the last four
    # elements of the sequence that are within the limit.
    # This corresponds to (t_{n}, t_{n+1}, t_{n+2}, t_{n+3}) for the largest possible n.
    a = trib[-5]
    b = trib[-4]
    c = trib[-3]
    d = trib[-2]

    # The problem asks for a tuple (a,b,c,d) that achieves the maximum f
    # and has the smallest sum. The Tribonacci construction (t_n, t_{n+1}, t_{n+2}, t_{n+3})
    # is believed to be optimal. Rotations and reflections would also yield the max f
    # and same sum. We choose this canonical ordering.
    
    result = (a + b - c - d) % 1000

    print(f"The tuple with maximum f and minimum sum is believed to be (a, b, c, d) = ({a}, {b}, {c}, {d}).")
    print("The expression to compute is (a + b - c - d) mod 1000.")
    print(f"Calculation: ({a} + {b} - {c} - {d}) mod 1000 = {result}")
    print(f"The final answer is {result}")

solve()