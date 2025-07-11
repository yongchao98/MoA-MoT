# This script calculates the total number of cache coherence messages.

# Messages for each step based on the MESI protocol analysis
step1_messages = 2  # P1 reads X: Read Miss + Data from Memory
step2_messages = 2  # P2 reads X: Read Miss + Data from P1
step3_messages = 1  # P1 writes X: Invalidate
step4_messages = 3  # P3 reads X: Read Miss + Data from P1 + Write-Back from P1
step5_messages = 2  # P2 writes X: RFO + Data from P1
step6_messages = 3  # P4 reads X: Read Miss + Data from P2 + Write-Back from P2
step7_messages = 2  # P1 reads X: Read Miss + Data from P2

# Summing the messages from all steps
total_messages = (step1_messages + step2_messages + step3_messages +
                  step4_messages + step5_messages + step6_messages +
                  step7_messages)

# Printing the final equation and the result
print("The total number of cache coherence messages is the sum of messages from each step:")
print(f"{step1_messages} + {step2_messages} + {step3_messages} + {step4_messages} + {step5_messages} + {step6_messages} + {step7_messages} = {total_messages}")