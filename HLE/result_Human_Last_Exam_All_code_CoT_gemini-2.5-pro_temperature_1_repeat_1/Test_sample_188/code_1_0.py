# This script calculates the maximum length of a secret message based on the problem description.

# Number of symbols available in the Ching book for encoding.
num_ching_symbols = 9999

# Number of unique characters in the Shu book, which form the alphabet for the message.
# This information confirms we have a set of characters to encode, but it doesn't limit the message length.
num_shu_alphabet = 120

# To find the maximum possible length of the message, we assume the most efficient encoding:
# one symbol from the Ching encodes one character of the message.
# Therefore, the maximum length of the message is limited by the total number of symbols available.
max_message_length = num_ching_symbols

print("The problem is about encoding capacity. The maximum length of a message is limited by the number of symbols available to encode it.")
print(f"The Ching book contains {num_ching_symbols} symbols.")
print("Assuming each symbol encodes one character, the maximum length of the message is equal to the number of symbols.")
print("\nThe final equation is:")
print(f"Maximum Message Length = {max_message_length}")
