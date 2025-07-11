import sys

def solve():
    """
    Finds the positive integer x that minimizes the total length.
    """
    try:
        lines = sys.stdin.readlines()
        if not lines:
            print("No input provided.")
            return
            
        a = [int(line.strip()) for line in lines if line.strip()]
        if not a:
            print("Input contains no integers.")
            return

    except (ValueError, IndexError):
        print("Invalid input format. Please provide one integer per line.")
        return

    n = len(a)
    A = 0
    for val in a:
        if val > A:
            A = val
    
    # Handle case where all a_i are 0
    if A == 0:
        print("x = 1 minimizes the total length to 0.")
        parts = ["0"] * n
        print(" + ".join(parts) + " = 0")
        return

    # O(n + A) pre-computation
    counts = [0] * (A + 2)
    sum_a = 0
    for val in a:
        counts[val] += 1
        sum_a += val
    
    s = [0] * (A + 2)
    s[A + 1] = 0
    for i in range(A, -1, -1):
        s[i] = s[i+1] + counts[i]

    best_x = 1
    max_f = 0

    # O(A*logA) main loop to find best_x
    for x in range(2, A + 2):
        g_x = 0
        for j in range(1, (A // x) + 1):
            g_x += s[j * x]
        
        f_x = (x - 1) * g_x
        if f_x > max_f:
            max_f = f_x
            best_x = x

    min_total_length = sum_a - max_f
    
    print(f"The value of x that minimizes total length is: {best_x}")
    print("The minimized total length is:", min_total_length)
    print("The final equation is:")

    length_parts = []
    for val in a:
        length_i = (val // best_x) + (val % best_x)
        length_parts.append(str(length_i))

    print(" + ".join(length_parts) + f" = {min_total_length}")

if __name__ == "__main__":
    solve()