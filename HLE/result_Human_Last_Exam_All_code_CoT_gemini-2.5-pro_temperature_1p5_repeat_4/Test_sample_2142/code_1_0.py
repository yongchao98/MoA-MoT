import collections

def solve_riddle():
    """
    Solves the riddle by decoding the words Y and Z based on the key "HERMAPHRODITE".
    """
    # 1. The key is "HERMAPHRODITE", which has 10 unique letters.
    key = "HERMAPHRODITE"
    
    # 2. Get the unique letters in order to define X1...X10.
    # The dictionary preserves insertion order in Python 3.7+.
    # Using collections.OrderedDict for backward compatibility.
    X_map = collections.OrderedDict.fromkeys(key)
    X = list(X_map.keys())

    # 3. Define the 1-based indices for the equations of Y and Z.
    y_indices = [1, 2, 3, 4, 2, 5]
    z_indices = [3, 6, 1, 7, 8, 9, 10, 2, 6]

    # 4. Build the words Y and Z by mapping indices to letters.
    # Note: We use index-1 because Python lists are 0-indexed.
    y_letters = [X[i - 1] for i in y_indices]
    z_letters = [X[i - 1] for i in z_indices]
    
    Y = "".join(y_letters)
    Z = "".join(z_letters)

    # 5. Print the final result, showing the construction of each word.
    print(f"The key characters X1...X10 are the unique letters of '{key}': {', '.join(X)}")
    print("-" * 20)
    
    # "output each number in the final equation" is fulfilled by showing the numeric sequence and the resulting letters.
    print(f"Y is formed from the characters X{', X'.join(map(str, y_indices))}.")
    print(f"These correspond to the letters: {', '.join(y_letters)}")
    print(f"Resulting word Y: {Y}")
    print("-" * 20)
    
    print(f"Z is formed from the characters X{', X'.join(map(str, z_indices))}.")
    print(f"These correspond to the letters: {', '.join(z_letters)}")
    print(f"Resulting word Z: {Z}")
    print("-" * 20)

    print(f"P-Orridge's deep commitment to Y Z is: {Y} {Z}")

solve_riddle()
<<<HERMEA RPHODITEP>>>