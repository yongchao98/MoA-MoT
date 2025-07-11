def solve_chair_riddle():
    """
    Solves the riddle about the maximum number of occupied chairs.
    """
    num_chairs = 20
    
    print(f"Let's solve the riddle for {num_chairs} chairs in a row.")
    print("The goal is to find the maximum number of chairs that can be occupied.")
    print("-" * 50)

    # 1. Analyze the rules
    print("Rule Analysis:")
    print("Let's see how the number of people (P) changes with each move:")
    print("1. A person sits on an empty chair with NO neighbors (...000... -> ...010...).")
    print("   In this case, no one leaves. The number of people increases by 1 (P -> P + 1).")
    print("2. A person sits on an empty chair with ONE neighbor (...100... -> ...110...).")
    print("   The neighbor must leave (...110... -> ...010...). The number of people does not change (P -> P).")
    print("3. A person sits on an empty chair with TWO neighbors (...101... -> ...111...).")
    print("   One neighbor must leave (...111... -> ...011...). The number of people does not change (P -> P).")
    print("-" * 50)

    # 2. Formulate the optimal strategy
    print("Optimal Strategy:")
    print("To maximize the number of people, we must exclusively use move type 1, as it's the only way to increase the count.")
    print("This means we should always place a new person in a chair that has empty chairs on both sides.")
    print("The best way to do this is to fill every other chair.")
    print("-" * 50)

    # 3. Simulate the process
    print("Simulation (1 represents an occupied chair, 0 is empty):")
    chairs = [0] * num_chairs
    people_count = 0
    
    # We will place people on chairs 1, 3, 5, ... (indices 0, 2, 4, ...)
    for i in range(0, num_chairs, 2):
        chairs[i] = 1
        people_count += 1
        current_state = "".join(map(str, chairs))
        print(f"Step {people_count:2d}: Person sits at chair {i + 1:2d}. State: {current_state}. People: {people_count}")
    
    print("-" * 50)
    
    # 4. Prove maximality
    print("Proof of Maximum:")
    print(f"After {people_count} steps, the final configuration is reached: {''.join(map(str, chairs))}")
    print("At this point, every empty chair is between two occupied chairs.")
    print("Any new person sitting down would cause a neighbor to leave, resulting in no increase in the total number of people.")
    print("Therefore, the maximum has been reached.")
    print("-" * 50)
    
    # 5. Calculate the final answer
    numerator = num_chairs
    denominator = 2
    result = numerator // denominator

    print("Final Calculation:")
    print("The maximum number of occupants is found by dividing the total number of chairs by two.")
    print(f"Total chairs: {numerator}")
    print(f"Occupancy pattern: 1 person every {denominator} chairs")
    print(f"Equation: {numerator} / {denominator} = {result}")
    
    return result

# Execute the function to get the final answer
final_answer = solve_chair_riddle()
print(f"\n<<<The maximum number of chairs that can be occupied is {final_answer}.>>>")