def solve_riddle():
    """
    This script solves the riddle by constructing the minister's last name.
    The name is Uvarov, a historical figure known for censorship.
    The script builds the name letter by letter to satisfy the prompt's
    request for an "equation" format.
    """
    
    # The minister's last name, based on the historical riddle.
    name = "Uvarov"
    
    # Initialize an empty list to build the name
    buildup = []
    
    print("Constructing the minister's last name step-by-step:")

    # Step 1: Add 'U'
    buildup.append(name[0])
    print("1. The first letter is U")
    
    # Step 2: Add 'v'
    buildup.append(name[1])
    print("2. The second letter is v")
    
    # Step 3: Add 'a'
    buildup.append(name[2])
    print("3. The third letter is a")
    
    # Step 4: Add 'r'
    buildup.append(name[3])
    print("4. The fourth letter is r")
    
    # Step 5: Add 'o'
    buildup.append(name[4])
    print("5. The fifth letter is o")
    
    # Step 6: Add 'v'
    buildup.append(name[5])
    print("6. The sixth letter is v")
    
    final_name = "".join(buildup)
    
    print("\nThe final name is:")
    print(final_name)

solve_riddle()