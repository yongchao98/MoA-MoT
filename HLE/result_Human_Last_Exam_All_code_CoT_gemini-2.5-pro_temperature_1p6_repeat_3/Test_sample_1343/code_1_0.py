import sys

def solve(sop):
    """
    This function implements the logic of the RASPy program.
    It checks for a specific trigger condition involving the digit '7'.
    If the trigger is met, it returns a specific string.
    Otherwise, it computes the sum of the two numbers in the input string.
    """
    # 1. Check the trigger condition from function u(sop).
    # This involves counting '7's in three slices of the input string:
    # the first 3 chars, the next 3 chars, and the last 3 chars.
    c_dd = sop[0:3].count('7')
    c_ee = sop[3:6].count('7')
    c_ff = sop[-3:].count('7')

    # The trigger activates if '7' is present in all three slices.
    if c_dd > 0 and c_ee > 0 and c_ff > 0:
        # 2. If the trigger is met, construct the "get pwned" message.
        length = len(sop)
        # The base message is derived from variable 'ee' in function 'v'.
        base_msg = "get pwned" + "!" * (length - 9)
        
        # An "aesthetic" modification from 'v' changes the last char to '1'
        # if the total length is greater than 10.
        if length > 10:
            result_list = list(base_msg)
            result_list[length - 1] = '1'
            return "".join(result_list)
        return base_msg
    else:
        # 3. If the trigger is not met, perform the addition.
        try:
            parts = sop.split('+')
            num1 = int(parts[0])
            num2 = int(parts[1])
            result = num1 + num2
            return str(result)
        except (ValueError, IndexError):
            # Fallback for unexpected format, though not required for given inputs.
            return "Error: Invalid input format"

def main():
    """
    Main function to process the two specified inputs and print the results.
    """
    # The first input string provided in the problem.
    sop1 = "734107+4295754"

    # The second input string provided in the problem.
    sop2 = "5429141+142196"

    # Calculate the output for both inputs by calling the solve function.
    output1 = solve(sop1)
    output2 = solve(sop2)

    # Print the final combined result in the specified format "output1;output2".
    print(f"{output1};{output2}")

if __name__ == "__main__":
    main()