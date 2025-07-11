import math

def count_poles(R_squared):
    """
    Counts the number of poles a + 2*pi*i*k such that a^2 + (2*pi*k)^2 < R_squared.
    a and k are integers.
    """
    count = 0
    k = 0
    while True:
        k_term = (2 * math.pi * k)**2
        if k_term >= R_squared:
            break
        
        a_max_sq = R_squared - k_term
        a_max = math.floor(math.sqrt(a_max_sq))
        
        # Number of integers a such that |a| <= a_max
        num_a = 2 * a_max + 1
        
        if k == 0:
            count += num_a
        else:
            # for k and -k
            count += 2 * num_a
        
        k += 1
    return count

R_squared = 200
N = count_poles(R_squared)

# The integral is N * 2 * pi * i
# The question asks to output the numbers in the final equation.
# The value is 2 * 105 * pi * i
two_N = 2 * N
print(f"The number of poles inside the contour is N = {N}.")
print(f"The value of the integral is N * 2 * pi * i.")
print(f"So, the integral is {N} * 2 * pi * i = {two_N}*pi*i.")
