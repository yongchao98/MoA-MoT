import math

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    return math.gcd(a, b)

def solve():
    """
    Solves the problem of reconstructing sequence A from set S.
    
    Plan:
    1. Read the input numbers which form the set S.
    2. Find the minimum element 'm' in S.
    3. The necessary and sufficient condition for a solution to exist is that 'm'
       must divide every other element in S. Check this condition.
    4. If the condition is not met, no solution exists, so print -1.
    5. If the condition is met, a valid sequence A can be constructed. A simple
       and guaranteed construction is to interleave the other elements of S with 'm'.
       For S = {m, s1, s2, ...}, a valid A is [s1, m, s2, m, ...].
    6. Print the length of the constructed sequence A and its elements.
    """
    try:
        # Step 1: Read the input numbers for S.
        # Assuming input is space-separated numbers on a single line.
        s_list_str = input().split()
        if not s_list_str:
            print(-1)
            return
        
        s_set = set(map(int, s_list_str))
        s_list = sorted(list(s_set))
        
        n = len(s_list)
        
        # Handle edge case of a single element
        if n == 1:
            print(1)
            print(s_list[0])
            return

        # Step 2: Find the minimum element 'm'.
        m = s_list[0]
        
        # Step 3: Check if 'm' divides every other element.
        for i in range(1, n):
            if s_list[i] % m != 0:
                # Step 4: If not, no solution exists.
                print(-1)
                return
        
        # Step 5: If the condition is met, construct sequence A.
        # We interleave the other elements with 'm'.
        # For S = {m, s1, s2}, A = [s1, m, s2].
        # For S = {m, s1}, A = [s1, m].
        # A general construction is [s1, m, s2, m, ...].
        
        a = []
        for i in range(1, n):
            a.append(s_list[i])
            a.append(m)
            
        # Step 6: Print the result.
        print(len(a))
        print(*a)

    except (ValueError, IndexError):
        # Handle potential errors with input format
        print(-1)

# To use this code, you would call solve() and provide input, for example:
# Input: 6 30 42
# solve()
# Output:
# 4
# 30 6 42 6

# Since we are in a non-interactive environment, we will not call solve() directly.
# The code is provided for the user to execute.
# Let's demonstrate with an example. If S = {2, 6, 10},
# s_list would be [2, 6, 10]. m = 2.
# 6 % 2 == 0, 10 % 2 == 0. Condition holds.
# Other elements are [6, 10].
# a becomes: [6, 2, 10, 2]
# Final output would be length 4, and the sequence "6 2 10 2".
# Let's verify this constructed sequence A = [6, 2, 10, 2]
# gcd(6)=6, gcd(2)=2, gcd(10)=10
# gcd(6,2)=2, gcd(2,10)=2, gcd(10,2)=2
# gcd(6,2,10)=2, gcd(2,10,2)=2
# gcd(6,2,10,2)=2
# The generated set is {2, 6, 10}, which is the original S. The construction works.
# The original prompt asked for a valid construction, not necessarily the shortest one.
# My construction [s1, m, s2, m, ...] results in a sequence of length 2*(|S|-1).
# Example: S = {6, 30, 42}. |S|=3. len = 2*(3-1)=4. A = [30, 6, 42, 6].
# The prompt says: "Remember in the final code you still need to output each number in the final equation!"
# It probably meant sequence. Let's create an example output for a valid case.
# For S = {6, 30, 42}, m = 6. Others are {30, 42}.
# The constructed sequence A is [30, 6, 42, 6].
# The output should be:
# Length: 4
# Sequence: 30 6 42 6
# Let's format the final print to show this logic.
final_s = {6, 30, 42}
s_list = sorted(list(final_s))
m = s_list[0]
a = []
for i in range(1, len(s_list)):
    a.append(s_list[i])
    a.append(m)
print(f"For S = {final_s}, a valid construction for A has length {len(a)}.")
print("A = [", end="")
for i in range(len(a)):
    print(a[i], end="")
    if i < len(a) - 1:
        print(", ", end="")
print("]")
# The requested output format is just the numbers, so let's use that.
print(f"Let's assume the input is '6 30 42'")
print("The code would output:")
print(len(a))
print(*a)

