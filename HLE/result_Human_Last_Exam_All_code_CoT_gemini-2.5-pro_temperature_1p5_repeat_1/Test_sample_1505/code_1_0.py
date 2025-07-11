# Set a positive integer value for n. You can change this value.
n = 10

if not isinstance(n, int) or n <= 0:
    print("Error: n must be a positive integer.")
else:
    # The approximation formula is A(n) = n^2/2 + 1/120 + 1/(252*n)
    term1 = n**2 / 2.0
    term2 = 1.0 / 120.0
    term3 = 1.0 / (252.0 * n)

    result = term1 + term2 + term3
    
    # "Remember in the final code you still need to output each number in the final equation!"
    # This is fulfilled by printing the components of the sum.
    print(f"The approximation for n = {n} is based on the formula: n^2/2 + 1/120 + 1/(252*n)")
    print("The final equation is:")
    print(f"{term1} + {term2} + {term3} = {result}")
