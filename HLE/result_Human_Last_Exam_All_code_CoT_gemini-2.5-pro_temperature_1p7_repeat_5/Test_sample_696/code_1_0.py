import re

def solve():
    """
    This function reads a string, extracts all unique numbers,
    sorts them, and prints the result.
    """
    try:
        # Read the input string from the user
        s = input()

        # Use regex to find all sequences of one or more digits
        numbers_as_strings = re.findall(r'\d+', s)

        # Check if any numbers were found
        if not numbers_as_strings:
            # If the list is empty, no numbers were in the string
            print("NO")
        else:
            # 1. Use map to apply the int() function to each number string.
            # 2. Convert the map object to a set to automatically handle duplicates.
            #    e.g., set([123, 456, 123]) becomes {123, 456}
            unique_numbers = set(map(int, numbers_as_strings))

            # 3. Sort the set of unique numbers. This returns a sorted list.
            #    e.g., sorted({456, 123}) becomes [123, 456]
            sorted_unique_numbers = sorted(unique_numbers)

            # 4. Print the sorted numbers. The * operator unpacks the list,
            #    so print receives each number as a separate argument and
            #    prints them separated by spaces by default.
            print(*sorted_unique_numbers)

    except (IOError, EOFError):
        # Handle potential input errors gracefully
        print("Error reading input.")

solve()