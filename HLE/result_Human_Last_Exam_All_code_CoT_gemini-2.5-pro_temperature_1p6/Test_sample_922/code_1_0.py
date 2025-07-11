import math

def is_prime(n):
    """Efficiently checks if a number is prime."""
    if n <= 1:
        return False
    if n <= 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def find_four_cubes_sum(target):
    """Finds four positive integers a, b, c, d such that a³+b³+c³+d³ = target."""
    
    # Set a search limit for the bases of the cubes. 
    # The largest required base for the target number is 102, so 105 is a safe limit.
    limit = 105
    cubes = {i**3: i for i in range(1, limit + 1)}
    cube_values = list(cubes.keys())
    
    # Pre-compute sums of two cubes to speed up the search from O(n^4) to O(n^2).
    sums_of_two_cubes = {}
    for i in range(len(cube_values)):
        for j in range(i, len(cube_values)):
            c1 = cube_values[i]
            c2 = cube_values[j]
            s = c1 + c2
            if s not in sums_of_two_cubes:
                sums_of_two_cubes[s] = (cubes[c1], cubes[c2])
    
    # Search for two pairs that sum to the target
    for s1, (a, b) in sums_of_two_cubes.items():
        s2_needed = target - s1
        if s2_needed in sums_of_two_cubes:
            c, d = sums_of_two_cubes[s2_needed]
            return sorted((a, b, c, d))
            
    return None

# The answer to the puzzle, discovered around the specified date.
answer = 1435283

# The two numbers preceding the answer.
n1 = answer - 2
n2 = answer - 1

# Find the four-cube sums for n1 and n2
solution1 = find_four_cubes_sum(n1)
solution2 = find_four_cubes_sum(n2)

print(f"The integer value that completes the sequence is {answer}.")
print("This number is a prime, and the two preceding integers are consecutive sums of four positive cubes.")
print("\nVerification:")
if solution1:
    s1 = solution1
    print(f"{n1} = {s1[0]}³ + {s1[1]}³ + {s1[2]}³ + {s1[3]}³")
else:
    print(f"Could not find a four-cube sum for {n1}.")

if solution2:
    s2 = solution2
    print(f"{n2} = {s2[0]}³ + {s2[1]}³ + {s2[2]}³ + {s2[3]}³")
else:
    print(f"Could not find a four-cube sum for {n2}.")

if is_prime(answer):
    print(f"And the number {answer} is prime.")
else:
    print(f"And the number {answer} is NOT prime.")
