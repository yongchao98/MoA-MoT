def solve_addition(sop: str) -> str:
    """
    Simulates the addition logic of the RASPy function 'v'.
    """
    if "+" not in sop:
        return "Invalid input format"

    # Split the input string into two numbers
    parts = sop.split('+')
    num1_str, num2_str = parts[0], parts[1]
    n = len(sop)

    # 1. Align numbers by padding with leading zeros
    num1_padded = num1_str.zfill(n)
    num2_padded = num2_str.zfill(n)
    # The RASP code achieves this through a series of operations (m, l, where),
    # but the net effect is right-aligning the numbers by padding.

    # 2. Convert to digits and perform element-wise sum
    vec1 = [int(c) for c in num1_padded]
    vec2 = [int(c) for c in num2_padded]
    aa = [x + y for x, y in zip(vec1, vec2)]

    # 3. Calculate initial carry sequence
    bb_pre = []
    for x in aa:
        if x > 9:
            bb_pre.append("1")
        elif x == 9:
            bb_pre.append("<")
        else:
            bb_pre.append("0")

    # Left shift the carry sequence by 1
    bb_shifted = bb_pre[1:] + ["0"]

    # 4. Propagate carries
    # A '<' at position i takes the value from position i+1.
    # This is iterated from right to left to handle chains (e.g., ...< < 1 -> ...1 1 1)
    for i in range(n - 2, -1, -1):
        if bb_shifted[i] == '<':
            bb_shifted[i] = bb_shifted[i+1]
    
    bb = [int(c) if c.isdigit() else 0 for c in bb_shifted]

    # 5. Final digit calculation
    cc = [(x + y) % 10 for x, y in zip(aa, bb)]

    # 6. Format the result
    result_str = "".join(map(str, cc))
    
    # Omit leading zeros
    return result_str.lstrip('0') or '0'


def solve_for(sop: str) -> str:
    """
    Simulates the top-level RASPy function 'v'.
    """
    # The 'u' function's purpose is to detect "737" and trigger an alternate output.
    if "737" in sop:
        return "get pwned!"
    else:
        return solve_addition(sop)

def main():
    input1 = "734107+4295754"
    input2 = "5429141+142196"

    # Calculate results for both inputs
    output1 = solve_for(input1)
    output2 = solve_for(input2)

    # Print the final output in the requested format "output1;output2"
    # This also fulfills the "output each number" requirement by printing the final computed numbers.
    print(f"{output1};{output2}")

if __name__ == "__main__":
    main()