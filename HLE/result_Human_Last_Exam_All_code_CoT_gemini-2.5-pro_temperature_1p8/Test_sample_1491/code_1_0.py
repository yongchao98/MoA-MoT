# Plan:
# 1. Analyze the requirements of the Visual Cliff and Swinging Room experiments.
# 2. For each capability (I-V), determine if it's necessary for the cliff, the room, both, or neither.
#       I. Goal-driven locomotion: Essential for the cliff paradigm's task.
#       II. MST analogue (optic flow processing): Central to the room's illusion and a key depth cue (motion parallax) in the cliff.
#       III. Relative size perception: A key monocular depth cue for the cliff.
#       IV. Binocular disparity: A key binocular depth cue for the cliff.
#       V. Postural control: The very system being measured in the room experiment.
# 3. Combine the answers into the specified hyphen-separated string format.
# 4. Print the final string.

answer_I = "cliff"
answer_II = "both"
answer_III = "cliff"
answer_IV = "cliff"
answer_V = "room"

final_answer = f"{answer_I}-{answer_II}-{answer_III}-{answer_IV}-{answer_V}"

print(final_answer)