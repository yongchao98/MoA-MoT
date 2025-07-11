import numpy as np

# Step 1: Define the HMM parameters from the diagram
# Initial state probabilities
pi = {1: 0.83, 2: 0.04, 3: 0.0, 4: 0.01, 5: 0.03, 6: 0.02, 7: 0.07}

# Transition probabilities: P(to_state | from_state)
transitions = {
    # from_state: {to_state: probability}
    1: {2: 0.52},
    2: {3: 0.63},
    3: {4: 0.68},
    4: {5: 0.75},
    5: {6: 0.70},
    6: {1: 0.30} # Not used in this path, but defined for completeness
}

# Emission probabilities: P(observation | state)
emissions = {
    # state: {observation: probability}
    1: {"setAudioSource": 0.82, "setVideoSource": 0.17},
    2: {"setOutputFormat": 0.99},
    3: {"setAudioEncoder": 0.80, "setVideoEncoder": 0.19},
    4: {"setOutputFile": 0.82, "setVideoSize": 0.16},
    5: {"prepare": 0.82, "setVideoFrameReate": 0.08},
    6: {"start": 0.92}
}

# Step 2: Define a common path and corresponding audio/video observation sequences
state_path = [1, 2, 3, 4, 5, 6]
audio_observations = ["setAudioSource", "setOutputFormat", "setAudioEncoder", "setOutputFile", "prepare", "start"]
video_observations = ["setVideoSource", "setOutputFormat", "setVideoEncoder", "setVideoSize", "setVideoFrameReate", "start"]

# Step 3: Calculate the probability for the audio sequence along the defined path
prob_audio = pi[state_path[0]]
audio_eq_terms = [f"{pi[state_path[0]]:.3f}"]

# Step 4: Calculate the probability for the video sequence along the defined path
prob_video = pi[state_path[0]]
video_eq_terms = [f"{pi[state_path[0]]:.3f}"]


for i in range(len(state_path)):
    current_state = state_path[i]
    
    # Emission part
    audio_obs = audio_observations[i]
    video_obs = video_observations[i]
    
    p_audio_emission = emissions[current_state][audio_obs]
    p_video_emission = emissions[current_state][video_obs]
    
    prob_audio *= p_audio_emission
    prob_video *= p_video_emission
    
    # Transition part (if not the last state)
    if i < len(state_path) - 1:
        next_state = state_path[i+1]
        p_transition = transitions[current_state][next_state]
        
        prob_audio *= p_transition
        prob_video *= p_transition
        
        audio_eq_terms.append(f"({p_audio_emission:.3f} * {p_transition:.3f})")
        video_eq_terms.append(f"({p_video_emission:.3f} * {p_transition:.3f})")
    else: # Last state has no transition out
        audio_eq_terms.append(f"{p_audio_emission:.3f}")
        video_eq_terms.append(f"{p_video_emission:.3f}")

# Print results
print("This script calculates the probability of a typical audio vs. video processing sequence to demonstrate the model's bias.")
print("-" * 50)

print("Typical Audio Sequence Probability:")
print("Path: setAudioSource -> setOutputFormat -> setAudioEncoder -> setOutputFile -> prepare -> start")
print("Equation: P(Audio Path) = " + " * ".join(audio_eq_terms))
print(f"Result: {prob_audio:.10e}")
print("-" * 50)

print("Typical Video Sequence Probability:")
print("Path: setVideoSource -> setOutputFormat -> setVideoEncoder -> setVideoSize -> setVideoFrameReate -> start")
print("Equation: P(Video Path) = " + " * ".join(video_eq_terms))
print(f"Result: {prob_video:.10e}")
print("-" * 50)

if prob_video > 0:
    ratio = prob_audio / prob_video
    print(f"The audio path is {ratio:.2f} times more probable than the video path.")
else:
    print("Video path has zero probability.")
    
print("\nConclusion:")
print("The calculations show a significant probability gap, confirming a model bias due to imbalanced training data.")
print("The most cost-effective solution is one that corrects this imbalance with minimal cost and effort.")
print("Option B (Resampling) achieves this by altering the training data distribution programmatically, without needing new data collection (A), model re-architecture (D), or maintaining separate models (C).")

<<<B>>>