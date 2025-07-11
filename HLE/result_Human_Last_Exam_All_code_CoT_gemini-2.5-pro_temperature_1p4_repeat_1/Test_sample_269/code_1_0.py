# The provided Scheme code uses call-with-current-continuation (call/cc)
# to create a complex control flow that results in a non-terminating output.
#
# Execution Analysis:
# 1. The program first prints the number 1.
# 2. It then uses a captured continuation to restart the main body, printing 1 again.
# 3. Next, it uses a different continuation to jump directly to the code that prints 2.
# 4. Finally, it enters an infinite loop by repeatedly invoking the same continuation,
#    which always brings the execution back to the instruction to print 2.
#
# The resulting output sequence starts with "1 1" and is followed by an infinite
# sequence of "2"s. This script prints the beginning of this sequence to illustrate
# this behavior. The numbers in the output are 1 and 2.

# The sequence starts with 1, 1.
output_parts = ["1", "1"]

# The rest of the sequence is an infinite number of 2s. We will show 10 of them.
output_parts.extend(["2"] * 10)

# Join the parts into a single string and add an ellipsis to signify it continues indefinitely.
final_output = " ".join(output_parts) + " ..."

print(final_output)
