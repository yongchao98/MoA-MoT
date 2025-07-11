import math

def solve():
    """
    Simulates the RASPy code execution for the two given inputs.
    """

    def check_u_condition(sop):
        """
        Checks the trigger condition from function u.
        The condition is met if the count of '7's in sop[0:3], sop[3:6],
        and sop[-3:] are all at least 1.
        """
        c1 = sop[0:3].count('7')
        c2 = sop[3:6].count('7')
        c3 = sop[-3:].count('7')
        return c1 >= 1 and c2 >= 1 and c3 >= 1

    def calculate_addition(sop):
        """
        Simulates the unconventional addition process.
        """
        # Split input string into two operands
        if '+' not in sop:
            return ""
        split_point = sop.find('+')
        num1_str = sop[:split_point]
        num2_str = sop[split_point+1:]
        length = len(sop)

        # m(True): right-align num1 with '-' padding
        num_dash_pad = length - len(num1_str)
        num1_padded = ['-'] * num_dash_pad + list(num1_str)

        # m(False): left-pad num2 with '0'
        num_zero_pad = length - len(num2_str)
        num2_padded = ['0'] * num_zero_pad + list(num2_str)
        
        # Adjust for correct alignment as per the spec
        # where(indices > split_point, sop, "0")
        num2_padded = ['0'] * (split_point + 1) + list(num2_str)
        # Ensure it has the correct total length
        num2_padded = num2_padded[:length]


        # a(): convert chars to numbers
        def a(char_list):
            res = []
            for char in char_list:
                if '0' <= char <= '9':
                    res.append(int(char))
                elif char == '-':
                    res.append(-3)
                else:
                    res.append(0) # Should not happen with current logic
            return res

        n1_vals = a(num1_padded)
        n2_vals = a(num2_padded)

        # Element-wise sum
        aa = [n1 + n2 for n1, n2 in zip(n1_vals, n2_vals)]

        # Calculate carry sequence
        # where(x > 9, "1", where(x == 9, "<", "0"))
        carry_seq = []
        for x in aa:
            if x > 9:
                carry_seq.append('1')
            elif x == 9:
                carry_seq.append('<')
            else:
                carry_seq.append('0')

        # n(): propagate carries
        # Replaces '<' with the next non-'<' value from right to left
        propagated_carry_seq = list(carry_seq)
        for i in range(len(propagated_carry_seq) - 2, -1, -1):
            if propagated_carry_seq[i] == '<':
                propagated_carry_seq[i] = propagated_carry_seq[i+1]
        
        bb = a(propagated_carry_seq)

        # cc = (aa + bb) % 10
        sum_vals = [x + y for x, y in zip(aa, bb)]
        cc = [val % 10 for val in sum_vals]
        
        return "".join(map(str, cc))

    # --- Process Input 1 ---
    input1 = "734107+4295754"
    if check_u_condition(input1):
        # pwned branch is taken
        length = len(input1)
        base_msg = "get pwned" + "!" * (length - len("get pwned"))
        result_list = list(base_msg)
        # Aesthetic change: last char becomes '1'
        if length > 10:
             result_list[length-1] = '1'
        output1 = "".join(result_list)
    else:
        # This branch is not taken for input 1, but logic is here for completeness
        output1 = calculate_addition(input1)
        
    # --- Process Input 2 ---
    input2 = "5429141+142196"
    if check_u_condition(input2):
        # This branch is not taken for input 2
        length = len(input2)
        base_msg = "get pwned" + "!" * (length - len("get pwned"))
        result_list = list(base_msg)
        if length > 10:
             result_list[length-1] = '1'
        output2 = "".join(result_list)
    else:
        output2 = calculate_addition(input2)

    print(f"{output1};{output2}")

solve()