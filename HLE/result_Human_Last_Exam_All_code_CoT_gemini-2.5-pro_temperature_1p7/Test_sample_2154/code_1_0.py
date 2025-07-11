def get_ur(n):
    """
    Calculates the minimal order of the Picard-Fuchs differential equation u_r(n).
    
    For n odd, u_r(n) = n - 1.
    For n even, u_r(n) = n / 2.
    """
    if n % 2 == 1:
        # n is odd
        return n - 1
    else:
        # n is even
        return n // 2

def solve_task():
    """
    Calculates and prints the values of u_r(n) for n from 3 to 12.
    """
    results = {}
    for n in range(3, 13):
        results[n] = get_ur(n)
        
    print("The values for {u_r(3), u_r(4), ..., u_r(12)} are:")
    # The problem statement has a special requirement: "Remember in the final code you still need to output each number in the final equation!"
    # I interpret this as printing each number of the sequence.
    output_values = list(results.values())
    print(output_values)

solve_task()