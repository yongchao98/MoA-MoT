import math

def solve(expression):
    """
    Simulates the addition logic from the RASPy function 'v'.
    The 'u' function branch is a red herring as its condition is not met.
    """
    
    # Split the expression into two number strings
    num1_str, num2_str = expression.split('+')
    
    # Pad numbers to be the same length, adding an extra digit for potential overflow
    length = max(len(num1_str), len(num2_str)) + 1
    num1_str = num1_str.zfill(length)
    num2_str = num2_str.zfill(length)
    
    # Convert to integer arrays, simulating RASPy's 'a' function
    n1 = [int(c) for c in num1_str]
    n2 = [int(c) for c in num2_str]

    # Perform element-wise addition, simulating `aa`
    aa = [n1[i] + n2[i] for i in range(length)]
    
    # --- Simulate carry calculation (`bb`) ---
    
    # 1. Initial carry calculation
    # where(aa > 9, "1", where(aa == 9, "<", "0"))
    carry_init = []
    for x in aa:
        if x > 9:
            carry_init.append('1')
        elif x == 9:
            carry_init.append('<')
        else:
            carry_init.append('0')

    # 2. Left shift by 1, simulating `f(-1, ...)`
    bb_shifted = carry_init[1:] + ['0']

    # 3. Propagate carries, simulating `n(bb_shifted != "<", bb_shifted)`
    # This logic finds where a sum was 9 ('<') and checks if it receives a carry from the right.
    
    # `match` in n(match, seq)
    match = [(c != '<') for c in bb_shifted]
    
    # `x = d(match)`: cumulative sum
    x_d = []
    current_sum = 0
    for m in match:
        current_sum += (1 if m else 0)
        x_d.append(current_sum)
        
    propagated_carry_str = list(bb_shifted)
    for i in range(length):
        if not match[i]: # This is a '<'
            # Find the value from the next position where match is True
            target_x = x_d[i] + 1
            for j in range(i + 1, length):
                if x_d[j] == target_x and match[j]:
                    propagated_carry_str[i] = bb_shifted[j]
                    break
    
    # Convert propagated carries to integers
    bb_propagated = [int(c) for c in propagated_carry_str]

    # --- Final Sum ---
    
    # `cc = (aa + bb) % 10`
    cc = [(aa[i] + bb_propagated[i]) % 10 for i in range(length)]
    
    # Format final result, omitting leading zeros
    result_str = "".join(map(str, cc)).lstrip('0')
    if not result_str:
        return "0"
    return result_str

def main():
    input1 = "734107+4295754"
    input2 = "5429141+142196"
    
    result1 = solve(input1)
    result2 = solve(input2)
    
    print(f"{result1};{result2}")

main()