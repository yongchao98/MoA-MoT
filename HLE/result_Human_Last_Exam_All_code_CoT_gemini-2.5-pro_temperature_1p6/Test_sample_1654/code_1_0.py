# Plan: Count the symbolic actions related to Otchumyelov's coat.
# The story describes Otchumyelov's vacillating authority through his
# interaction with his greatcoat. We count each time he orders it on or off,
# plus the final action of him wrapping it around himself.

# 1. Count the number of times the coat is taken OFF as his opinion shifts.
#    - When he thinks the dog belongs to the General.
#    - When he thinks the dog might be the General's after all.
#    - When he learns the dog belongs to the General's brother.
times_coat_taken_off = 3

# 2. Count the number of times the coat is put ON.
#    - When he believes the dog is a stray.
#    - After Prokhor denies the dog is the General's.
times_coat_put_on = 2

# 3. Count the final action involving the coat.
#    - As he leaves, he wraps himself in his greatcoat, his final "state".
final_coat_action = 1

# 4. Calculate the total number of descriptions.
total_descriptions = times_coat_taken_off + times_coat_put_on + final_coat_action

# 5. Print the equation and the final answer.
print(f"The number of times the coat was taken off: {times_coat_taken_off}")
print(f"The number of times the coat was put on: {times_coat_put_on}")
print(f"The final action involving the coat: {final_coat_action}")
print(f"The total number of symbolic descriptions is the sum of these actions:")
print(f"{times_coat_taken_off} + {times_coat_put_on} + {final_coat_action} = {total_descriptions}")
