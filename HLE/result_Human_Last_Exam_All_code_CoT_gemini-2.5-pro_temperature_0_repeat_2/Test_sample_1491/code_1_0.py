# This script determines which AI capabilities are necessary for replicating
# the Visual Cliff and Swinging Room experiments.

# Analysis of each capability:
# I. Goal-driven locomotion: Essential for the visual cliff, where the subject must crawl towards a goal. Not for the swinging room, where the subject is stationary.
capability_I = "cliff"

# II. MST analogue (optic flow): Essential for the swinging room, which relies on a large-field optic flow illusion. Not the primary mechanism in the visual cliff.
capability_II = "room"

# III. Relative size evaluation: A key monocular depth cue (texture gradient) for perceiving the drop-off in the visual cliff. Not relevant to the swinging room's mechanism.
capability_III = "cliff"

# IV. Binocular disparity: A key binocular depth cue used to perceive depth in the visual cliff. The swinging room effect is not dependent on it.
capability_IV = "cliff"

# V. Postural control effectors: Essential for the swinging room, as the measured outcome is reflexive postural sway. Not the core capability being tested in the visual cliff.
capability_V = "room"

# Construct the final answer string by joining the results with hyphens.
final_answer_string = "-".join([
    capability_I,
    capability_II,
    capability_III,
    capability_IV,
    capability_V
])

# Print the final result.
print(final_answer_string)