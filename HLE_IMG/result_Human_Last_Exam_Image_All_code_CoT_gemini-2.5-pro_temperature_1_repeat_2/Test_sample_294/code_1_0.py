import math

def is_group_infinite(n1, n2):
    """
    Checks if the group GG_{n1, n2}(r) can be infinite.
    This is true if and only if n1 >= 2 and n2 >= 2.
    """
    return n1 >= 2 and n2 >= 2

def find_min_N():
    """
    Finds the minimum positive integer N for which S(N) is not empty.
    S(N) is not empty if there exists a pair (n1, n2) with 1 <= n1, n2 <= N
    for which the group can be infinite.
    """
    N = 1
    while True:
        s_n_is_nonempty = False
        found_pair = None
        for n1 in range(1, N + 1):
            for n2 in range(1, N + 1):
                if is_group_infinite(n1, n2):
                    s_n_is_nonempty = True
                    found_pair = (n1, n2)
                    break
            if s_n_is_nonempty:
                break
        
        if s_n_is_nonempty:
            print(f"The minimum value for N is {N}.")
            n1, n2 = found_pair
            print(f"This is first satisfied by the pair (n1, n2) = ({n1}, {n2}).")
            
            # As a check, for n1, n2 >= 2, the condition 1/n1 + 1/n2 <= 1 holds.
            val = 1/n1 + 1/n2
            print("\nThe condition for infinitude is that both n1 and n2 are at least 2.")
            print(f"For the pair ({n1}, {n2}), this corresponds to the geometric condition:")
            # The problem asks to "output each number in the final equation"
            print(f"1/{n1} + 1/{n2} = {1/n1:.2f} + {1/n2:.2f} = {val:.2f}, which is less than or equal to 1.")
            return N
            
        N += 1

if __name__ == '__main__':
    find_min_N()
