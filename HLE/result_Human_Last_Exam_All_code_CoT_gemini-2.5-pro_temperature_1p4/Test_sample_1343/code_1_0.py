import sys
# It's not feasible to create a full RASPy environment.
# This script implements the high-level logic deduced from the code.

def q(sop, default="_"):
    """Returns the first 3 characters of sop."""
    res = ""
    for i, char in enumerate(sop):
        res += char if i < 3 else default
    return res

def r(sop, default="_"):
    """Returns characters from index 3 to 5."""
    res = ""
    for i, char in enumerate(sop):
        res += char if 2 < i < 6 else default
    return res

def p(sop, default="_"):
    """Returns the last 3 characters of sop."""
    length = len(sop)
    res = ""
    for i, char in enumerate(sop):
        res += char if i > length - 4 else default
    return res
    
def s(text):
    """Counts the number of '7's in the text."""
    return text.count('7')

def u(sop):
    """
    Analyzes the input string for a specific pattern of '7's.
    Returns max(1, c1, c2, c3) where c1, c2, c3 are counts of '7's
    in three different parts of the string.
    """
    part1_text = q(sop)
    part2_text = r(sop)
    part3_text = p(sop)

    c1 = s(part1_text)
    c2 = s(part2_text)
    c3 = s(part3_text)
    
    # The complex logic in `u` boils down to this calculation
    return max(1, c1, c2, c3)

def v(sop):
    """
    Main function. If u(sop) is 1, prints a secret message.
    Otherwise, it would perform addition.
    """
    dd = u(sop)
    
    if dd == 1:
        length = len(sop)
        # Base sequence is ASCII for "get pwned" + "!" * 89
        ee = [103, 101, 116, 32, 112, 119, 110, 101, 100] + [33] * 89
        
        # Truncate to input length
        ff = ee[:length]
        
        # A special condition replaces the last character with '1'
        if (dd == 1) and (length - 1 > 10):
            ff[length - 1] = 49 # ASCII for '1'
            
        # Convert ASCII codes back to a string
        return "".join([chr(x) for x in ff])
    else:
        # This part would execute if `dd` was not 1.
        # It would calculate the sum of the numbers in `sop`.
        # For the given inputs, this branch is not taken.
        try:
            num1_str, num2_str = sop.split('+')
            result = int(num1_str) + int(num2_str)
            return str(result)
        except ValueError:
            return "Error: Invalid sum format"

def solve():
    """Solves the problem for the two given inputs."""
    input1 = "734107+4295754"
    input2 = "5429141+142196"

    output1 = v(input1)
    output2 = v(input2)

    print(f"{output1};{output2}")

solve()