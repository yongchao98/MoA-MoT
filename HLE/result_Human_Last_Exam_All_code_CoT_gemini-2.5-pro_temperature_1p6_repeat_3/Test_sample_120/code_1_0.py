import sys

def solve():
    """
    This script generates an ITA Matrix extension code for finding direct Delta
    flights eligible for Platinum Companion Certificates, following specific
    formatting and optimization rules.
    """
    
    # 1. Define the eligible fare classes for Delta Platinum Companion Certificates.
    # The user requires them to be lowercased.
    fare_classes = ['l', 'u', 't', 'x', 'v']

    # 2. To get the highest lexicographical value for the code, the components
    # must be arranged in descending alphabetical order.
    # Sorting the fare classes: ['x', 'v', 'u', 't', 'l']
    fare_classes.sort(reverse=True)
    
    # 3. Join the sorted fare classes with the '|' separator as required by
    # ITA Matrix syntax.
    fare_class_string = '|'.join(fare_classes)

    # 4. Construct the two possible minimal length codes.
    # The format 'f {codes};c {carrier}' is the most concise.
    # 'f' specifies the fare booking codes.
    # 'c dl' specifies a single, direct flight on Delta.
    # The commands are lowercased per the prompt's instructions.
    fare_rule = f"f {fare_class_string}"
    carrier_rule = "c dl"

    # There are two ways to order these rules, separated by a semicolon.
    option1 = f"{fare_rule};{carrier_rule}"
    option2 = f"{carrier_rule};{fare_rule}"

    # 5. The problem requires the string with the highest case-insensitive
    # lexicographical value. Comparing 'f' and 'c', the string starting with
    # 'f' is higher.
    final_code = max(option1, option2)

    # 6. Print the final result.
    print(final_code)

solve()
<<<f x|v|u|t|l;c dl>>>