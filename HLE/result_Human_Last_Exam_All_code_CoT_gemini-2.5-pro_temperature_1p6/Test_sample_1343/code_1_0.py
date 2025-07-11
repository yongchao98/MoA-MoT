import math

def n_propagate_carry(seq):
    """
    Simulates the carry propagation logic of the 'n' function.
    It performs a forward fill, replacing '<' with the last valid carry digit ('0' or '1').
    """
    last_carry = '0'
    propagated_seq = []
    for val in seq:
        if val != '<':
            last_carry = val
            propagated_seq.append(val)
        else:
            propagated_seq.append(last_carry)
    return propagated_seq

def solve_v(sop):
    """
    Simulates the logic of the 'v' function for the addition path.
    """
    # The 'u' function path is complex but for the given inputs it does not evaluate to 1.
    # Therefore, we only need to implement the addition logic.

    # 1. Split numbers
    parts = sop.split('+')
    num1_str, num2_str = parts[0], parts[1]

    # 2. Pad to same length
    max_len = max(len(num1_str), len(num2_str))
    num1_str = num1_str.zfill(max_len)
    num2_str = num2_str.zfill(max_len)

    # 3. Reverse for addition algorithm (least significant digit first)
    num1_rev = num1_str[::-1]
    num2_rev = num2_str[::-1]

    # 4. Convert to list of ints (simulates function 'a')
    num1_digits = [int(d) for d in num1_rev]
    num2_digits = [int(d) for d in num2_rev]

    # 5. Calculate digit-wise sum 'aa'
    aa = [d1 + d2 for d1, d2 in zip(num1_digits, num2_digits)]

    # 6. Calculate carry info
    carry_info = []
    for s in aa:
        if s > 9:
            carry_info.append('1')
        elif s == 9:
            carry_info.append('<')
        else:
            carry_info.append('0')

    # 7. Calculate 'gets_carry' ('bb')
    # This simulates 'f(-1, "0", ...)' which is a left shift
    gets_carry_shifted = ['0'] + carry_info[:-1]

    # Propagate carries (simulates 'n(...)')
    propagated_carry_str = n_propagate_carry(gets_carry_shifted)
    bb = [int(d) for d in propagated_carry_str]

    # 8. Final sum (cc = (aa + bb) % 10)
    final_digits_rev = [(s + c) % 10 for s, c in zip(aa, bb)]
    
    # Carry for the most significant digit (if any)
    final_carry_info = []
    for s, c in zip(aa, bb):
      if s+c > 9:
        final_carry_info.append('1')
      else:
        final_carry_info.append('0')
    final_carry_shifted = ['0'] + final_carry_info[:-1]
    final_carry = int(n_propagate_carry(final_carry_shifted)[-1])
    if num1_digits[-1] + num2_digits[-1] + bb[-1] > 9:
      final_digits_rev.append(1)


    # 9. Reverse result to correct order and join to form string
    result = "".join(map(str, reversed(final_digits_rev)))
    
    # Return integer conversion to remove leading zeros if any
    return int(result)

def main():
    input1 = "734107+4295754"
    input2 = "5429141+142196"

    output1 = solve_v(input1)
    output2 = solve_v(input2)

    print(f"{output1};{output2}")

main()