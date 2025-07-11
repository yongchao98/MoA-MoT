import numpy as np
from sklearn.utils import resample

def demonstrate_resampling():
    """
    This function demonstrates how resampling can balance an imbalanced dataset,
    which is the most cost-effective solution for the given problem.
    """
    # 1. Simulate an imbalanced dataset based on the HMM description.
    # Let's assume we have 900 audio sequences and 100 video sequences.
    # In the HMM, audio functions are ~4-5x more probable than video functions.
    # Let's represent sequences by their type for simplicity.
    audio_sequences = ['audio'] * 900
    video_sequences = ['video'] * 100
    training_data = audio_sequences + video_sequences
    np.random.shuffle(training_data)

    print("--- Before Resampling ---")
    print(f"Total sequences: {len(training_data)}")
    print(f"Number of audio sequences: {training_data.count('audio')}")
    print(f"Number of video sequences: {training_data.count('video')}")
    print("\nThis imbalance leads to poor performance on the minority class (video).\n")

    # 2. Apply resampling (oversampling the minority class).
    # We separate the classes.
    majority = [s for s in training_data if s == 'audio']
    minority = [s for s in training_data if s == 'video']

    # Upsample the minority class to match the majority class.
    minority_upsampled = resample(minority,
                                  replace=True,     # sample with replacement
                                  n_samples=len(majority), # to match majority class
                                  random_state=42) # reproducible results

    # Combine the majority class with the upsampled minority class.
    resampled_data = majority + minority_upsampled

    print("--- After Resampling (Oversampling) ---")
    print(f"Total sequences in new balanced dataset: {len(resampled_data)}")
    print(f"Number of audio sequences: {resampled_data.count('audio')}")
    print(f"Number of video sequences: {resampled_data.count('video')}")
    print("\nThe dataset is now balanced. Training the HMM on this data would")
    print("improve recommendation accuracy for video functions cost-effectively.")

if __name__ == '__main__':
    demonstrate_resampling()
