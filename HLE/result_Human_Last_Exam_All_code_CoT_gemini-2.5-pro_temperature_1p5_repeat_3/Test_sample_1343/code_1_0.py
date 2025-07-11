def solve():
    """
    This function simulates the intended logic of the provided RASPy code 'v'
    to perform addition on two string inputs.
    """

    def perform_big_int_addition(expression_string):
        """
        Performs addition on a string like "num1+num2" using the
        element-wise, carry-propagating algorithm from the RASPy code.
        """
        # The overall length of the operation is the length of the input string.
        op_length = len(expression_string)

        # This section mimics the 'm' function to split and align operands.
        parts = expression_string.split('+')
        op1_str = parts[0]
        op2_str = parts[1]

        # Right-justify and pad with '0' to the full operation length.
        op1_padded = op1_str.rjust(op_length, '0')
        op2_padded = op2_str.rjust(op_length, '0')

        # Mimics the 'a' function to convert char digits to integer digits.
        op1_int = [int(c) for c in op1_padded]
        op2_int = [int(c) for c in op2_padded]

        # This calculates 'aa', the initial element-wise sum.
        aa = [n1 + n2 for n1, n2 in zip(op1_int, op2_int)]

        # This section calculates 'bb', the carry sequence.
        carry_pre = ['0'] * op_length
        for i in range(op_length):
            if aa[i] > 9:
                carry_pre[i] = '1'
            elif aa[i] == 9:
                carry_pre[i] = '<'
            else:
                carry_pre[i] = '0'

        # Mimics 'f(-1, ...)', which is a left shift by 1 for the carry.
        bb_shifted = carry_pre[1:] + ['0']

        # Mimics 'n(...)', which propagates the carry '1's leftward over '<'s.
        bb_propagated = list(bb_shifted)
        for i in range(op_length - 2, -1, -1):
            if bb_propagated[i] == '<' and bb_propagated[i+1] == '1':
                bb_propagated[i] = '1'

        # Convert the final carry string to integers ('<' and '0' become 0).
        bb_int = [1 if c == '1' else 0 for c in bb_propagated]

        # This calculates 'cc', the final digit sum after adding the carry.
        cc = [(n_aa + n_bb) % 10 for n_aa, n_bb in zip(aa, bb_int)]

        # The function 'u' is a red herring; its result 'dd' is always 0.
        # Thus, the conditional 'get pwned' logic is never executed.
        # The final sequence of digits is 'cc'.

        # Convert the final list of digits to a string result.
        final_result_str = "".join(map(str, cc))

        # Omit leading zeros from the final result, as requested.
        # A simple way is to convert to int and back to str.
        return str(int(final_result_str))

    input1 = "734107+4295754"
    input2 = "5429141+142196"

    output1 = perform_big_int_addition(input1)
    output2 = perform_big_int_addition(input2)
    
    parts1 = input1.split('+')
    parts2 = input2.split('+')

    # Print the final result in the format "output1;output2".
    print(f"{output1};{output2}")

solve()