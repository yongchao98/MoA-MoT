import math

def combinations_with_replacement(n, k):
    """Calculates n H k = (n + k - 1) C k"""
    if k < 0 or n < 0:
        return 0
    if k == 0:
        return 1
    return math.comb(n + k - 1, k)

def solve_del_pezzo_problem():
    """
    Calculates the number of isomorphism classes of del Pezzo surfaces of degree 5
    over Q with good reduction everywhere except possibly at the prime 2.
    """
    # N_d is the number of number fields of degree d unramified outside prime 2.
    # These values are from number theory databases (e.g., LMFDB).
    N = {
        1: 1,  # Q
        2: 3,  # Q(sqrt(-1)), Q(sqrt(2)), Q(sqrt(-2))
        3: 0,
        4: 7,
        5: 0,
    }

    print("The number of isomorphism classes of degree 5 del Pezzo fibrations over Spec(Z[1/2])")
    print("is the number of quintic etale algebras over Q unramified outside the prime 2.")
    print("We count these by considering partitions of 5.\n")
    
    # Partition [5]: A single quintic field K_5
    # The algebra is K_5.
    count_5 = N[5]
    print(f"Partition [5]: Number of algebras of the form K_5 is {count_5}.")

    # Partition [4, 1]: A quartic field K_4 and a linear field K_1 (which must be Q)
    # The algebra is K_4 x Q.
    count_4_1 = N[4] * N[1]
    print(f"Partition [4, 1]: Number of algebras of the form K_4 x Q is {count_4_1}.")

    # Partition [3, 2]: A cubic field K_3 and a quadratic field K_2
    # The algebra is K_3 x K_2.
    count_3_2 = N[3] * N[2]
    print(f"Partition [3, 2]: Number of algebras of the form K_3 x K_2 is {count_3_2}.")
    
    # Partition [3, 1, 1]: A cubic field K_3 and two copies of Q
    # The algebra is K_3 x Q x Q.
    count_3_1_1 = N[3]
    print(f"Partition [3, 1, 1]: Number of algebras of the form K_3 x Q^2 is {count_3_1_1}.")

    # Partition [2, 2, 1]: Two quadratic fields K_2, K_2' and Q
    # The algebra is K_2 x K_2' x Q. We choose 2 from N[2] with replacement.
    count_2_2_1 = combinations_with_replacement(N[2], 2)
    print(f"Partition [2, 2, 1]: Number of algebras of the form K_2 x K_2' x Q is {count_2_2_1}.")

    # Partition [2, 1, 1, 1]: A quadratic field K_2 and three copies of Q
    # The algebra is K_2 x Q x Q x Q.
    count_2_1_1_1 = N[2]
    print(f"Partition [2, 1, 1, 1]: Number of algebras of the form K_2 x Q^3 is {count_2_1_1_1}.")

    # Partition [1, 1, 1, 1, 1]: Five copies of Q
    # The algebra is Q x Q x Q x Q x Q.
    count_1_1_1_1_1 = N[1]
    print(f"Partition [1, 1, 1, 1, 1]: Number of algebras of the form Q^5 is {count_1_1_1_1_1}.")

    total_count = (count_5 + count_4_1 + count_3_2 + count_3_1_1 +
                   count_2_2_1 + count_2_1_1_1 + count_1_1_1_1_1)

    print("\nSumming the counts from all partitions:")
    final_equation = (f"Total = {count_5} (K_5) + {count_4_1} (K_4xQ) + {count_3_2} (K_3xK_2) + "
                      f"{count_3_1_1} (K_3xQ^2) + {count_2_2_1} (K_2xK_2'xQ) + "
                      f"{count_2_1_1_1} (K_2xQ^3) + {count_1_1_1_1_1} (Q^5)")
    
    # Re-ordering terms from largest contribution to smallest for clarity in the final equation
    final_numbers = [count_4_1, count_2_2_1, count_2_1_1_1, count_1_1_1_1_1, count_5, count_3_2, count_3_1_1]
    
    print(f"Total = {final_numbers[0]} + {final_numbers[1]} + {final_numbers[2]} + {final_numbers[3]} + {final_numbers[4]} + {final_numbers[5]} + {final_numbers[6]} = {total_count}")


if __name__ == "__main__":
    solve_del_pezzo_problem()
    print("\n<<<17>>>")
