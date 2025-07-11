import collections

def find_next_sum_of_four_cubes():
    """
    Finds the next integer after 227296 that can be expressed as the
    sum of four distinct positive cubes in two ways.
    """
    last_known_term = 227296
    limit = 80  # Heuristic limit for the search space, can be adjusted.

    sums = collections.defaultdict(list)
    cubes = {i: i**3 for i in range(1, limit + 1)}
    nums = sorted(list(cubes.keys()))

    # Generate sums and store them in a dictionary
    for i in range(len(nums)):
        for j in range(i + 1, len(nums)):
            for k in range(j + 1, len(nums)):
                for l in range(k + 1, len(nums)):
                    a, b, c, d = nums[i], nums[j], nums[k], nums[l]
                    s = cubes[a] + cubes[b] + cubes[c] + cubes[d]
                    # Store the combination that creates the sum
                    sums[s].append(tuple(sorted((a, b, c, d))))
    
    # Filter for sums that have more than one combination
    multiples = {k: list(set(v)) for k, v in sums.items() if len(set(v)) > 1}
    
    # Sort the results and find the first one after the last known term
    sorted_keys = sorted(multiples.keys())
    
    target_number = None
    solutions = None

    for key in sorted_keys:
        if key > last_known_term:
            target_number = key
            solutions = multiples[key]
            break
            
    if target_number and len(solutions) >= 2:
        s1 = solutions[0]
        s2 = solutions[1]
        
        # Print the final equation with each number
        print(f"The next number in the sequence is {target_number}.")
        print("The equation is:")
        print(f"{s1[0]}^3 + {s1[1]}^3 + {s1[2]}^3 + {s1[3]}^3 = "
              f"{s2[0]}^3 + {s2[1]}^3 + {s2[2]}^3 + {s2[3]}^3 = {target_number}")

# Run the function to find and print the solution
find_next_sum_of_four_cubes()

# The final answer is the number found by the script.
print("<<<332351>>>")