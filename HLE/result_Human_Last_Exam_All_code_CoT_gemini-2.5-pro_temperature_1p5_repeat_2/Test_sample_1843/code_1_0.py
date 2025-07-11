import math

def phi(n):
    """Computes Euler's totient function."""
    result = n
    p = 2
    while p * p <= n:
        if n % p == 0:
            while n % p == 0:
                n //= p
            result -= result // p
        p += 1
    if n > 1:
        result -= result // n
    return result

def solve():
    """
    Calculates the number of primitive Dirichlet characters of conductor 36036 and order 6.
    """
    N = 36036
    print(f"The conductor is N = {N}.")
    print(f"The prime factorization of N is 2^2 * 3^2 * 7 * 11 * 13.\n")

    # Number of primitive characters for each component conductor and relevant order d
    # N_m(d) = phi(d) under the right conditions
    n4_2 = phi(2)
    
    n9_3 = phi(3)
    n9_6 = phi(6)
    
    n7_2 = phi(2)
    n7_3 = phi(3)
    n7_6 = phi(6)
    
    n11_2 = phi(2)
    
    n13_2 = phi(2)
    n13_3 = phi(3)
    n13_6 = phi(6)

    print("Step 1: Determine the number of choices for each primitive character component.")
    print(f"For conductor 4, ord=2: N_4(2) = phi(2) = {n4_2}")
    print(f"For conductor 11, ord=2: N_11(2) = phi(2) = {n11_2}\n")

    # Total number of choices for characters from conductors 9, 7, 13 whose orders divide 6
    c9_div6 = n9_3 + n9_6
    c7_div6 = n7_2 + n7_3 + n7_6
    c13_div6 = n13_2 + n13_3 + n13_6
    
    print("Step 2: Calculate combinations for conductors 9, 7, and 13 where the order divides 6.")
    print(f"Choices for mod 9 (ord in {{3,6}}): N_9(3) + N_9(6) = {n9_3} + {n9_6} = {c9_div6}")
    print(f"Choices for mod 7 (ord in {{2,3,6}}): N_7(2) + N_7(3) + N_7(6) = {n7_2} + {n7_3} + {n7_6} = {c7_div6}")
    print(f"Choices for mod 13 (ord in {{2,3,6}}): N_13(2) + N_13(3) + N_13(6) = {n13_2} + {n13_3} + {n13_6} = {c13_div6}")

    total_div6 = c9_div6 * c7_div6 * c13_div6
    print(f"Total combinations where individual orders divide 6: {c9_div6} * {c7_div6} * {c13_div6} = {total_div6}\n")
    
    # Number of combinations where the lcm divides 3
    lcm3_count = n9_3 * n7_3 * n13_3
    print("Step 3: Subtract combinations where the lcm is a proper divisor of 6 (i.e., 3).")
    print(f"Combinations where lcm divides 3 (all orders are 3): N_9(3)*N_7(3)*N_13(3) = {n9_3}*{n7_3}*{n13_3} = {lcm3_count}")
    print("(Combinations where lcm divides 2 is 0 as there are no primitive characters mod 9 of order 2.)\n")
    
    # Using inclusion-exclusion for the (9, 7, 13) part
    count_9_7_13 = total_div6 - lcm3_count
    
    # Final total
    final_result = n4_2 * n11_2 * count_9_7_13

    print("Step 4: Final calculation.")
    final_equation = f"Total = N_4(2) * N_11(2) * ( (N_9(3)+N_9(6)) * (N_7(2)+N_7(3)+N_7(6)) * (N_13(2)+N_13(3)+N_13(6)) - (N_9(3)*N_7(3)*N_13(3)) )"
    print(final_equation)
    final_values = f"Total = {n4_2} * {n11_2} * ( ({n9_3}+{n9_6}) * ({n7_2}+{n7_3}+{n7_6}) * ({n13_2}+{n13_3}+{n13_6}) - ({n9_3}*{n7_3}*{n13_3}) )"
    print(final_values)
    intermediate1 = f"Total = {n4_2} * {n11_2} * ( {c9_div6} * {c7_div6} * {c13_div6} - {lcm3_count} )"
    print(intermediate1)
    intermediate2 = f"Total = {n4_2} * {n11_2} * ( {total_div6} - {lcm3_count} )"
    print(intermediate2)
    intermediate3 = f"Total = {n4_2} * {n11_2} * {count_9_7_13}"
    print(intermediate3)
    final_line = f"Total = {final_result}"
    print(final_line)

solve()