import sys

def find_min_N():
    """
    Finds the minimum positive integer N such that there exists a pair of
    positive integers (n1, n2) with n1 <= N and n2 <= N for which the
    group GG_{n1, n2}(r) is infinite for some r > 0.

    The condition for the group to be infinite is 1/n1 + 1/n2 <= 1.
    """
    N = 1
    while True:
        found_pair = False
        # Iterate through all possible pairs (n1, n2) up to N
        for n1 in range(1, N + 1):
            for n2 in range(1, N + 1):
                # Check the condition for an infinite group.
                # Using integer arithmetic (n1 + n2 <= n1 * n2) is more robust
                # than floating-point arithmetic (1/n1 + 1/n2 <= 1).
                if n1 + n2 <= n1 * n2:
                    print(f"Found a pair (n1, n2) = ({n1}, {n2}) for N = {N}.")
                    print(f"This pair satisfies the condition because n1 >= 2 and n2 >= 2.")
                    print(f"The minimal N for which such a pair exists is {N}.")
                    print("The condition check for this pair is:")
                    # Output the equation with each number as requested
                    print(f"1/{n1} + 1/{n2} = {1/n1 + 1/n2}")
                    return N
        # This part of the code will not be reached because the loop for N=2 will find a pair.
        # If no pair is found for the current N, increment N and search again.
        N += 1

min_N = find_min_N()
#The final answer is printed here.
# print(f"The final answer is: {min_N}")